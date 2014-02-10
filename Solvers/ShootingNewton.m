classdef ShootingNewton < Solver
    %ShootingNewton.m Time domain periodic steady-state solver
    %
    % ShootingNewton properties:
    %   ShootingTolerance     - Tolerance of the periodic boundary condition
    % 	NewtonTolerance       - Tolerance of iteration occuring at each time point
    % 	GMRESTolerance        - Tolerance of the initial condition correction
    %  	MaxShootingIterations - Maximum number of outer iterations
    %  	MaxNewtonIterations   - Maximum number of iterations at each time point
    %  	MaxGMRESIterations    - Maximum number of iterations for the correction
    %  	RungeKuttaStages      - Number of stages used in the numerical integration
    % 	StoreDecompositions   - Toggle matrix decomposition storage for the GMRES stage
    %
    % ShootingNewton inherits properties and methods Solver.
    %
    % See also MotorProto, Solver,
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%NewtonTolerance - Sets the relative residual tolerance
 	%	NewtonTolerance sets the relative residual tolerance for the Newton-Raphson
    %   iteration that occurs at each time step. The default value is sqrt(eps).
    %
    % See also Static, MaxNewtonIterations
    NewtonTolerance;
    
 	%MaxNewtonIterations - Sets the maximum number of iterations
    %   MaxNewtonIterations sets the maximum number of Newton-Raphson iterations
    %   occuring at each time step. The default value is 100.
    %
    % See also Static, NewtonTolerancerations
    MaxNewtonIterations;
%}

    properties
        ShootingTolerance     = ShootingNewton.setProperty(1e-6);
        NewtonTolerance       = ShootingNewton.setProperty(1e-6);
        GMRESTolerance        = ShootingNewton.setProperty(1e-7);
        MaxShootingIterations = ShootingNewton.setProperty(100);
        MaxNewtonIterations   = ShootingNewton.setProperty(100);
        MaxGMRESIterations    = ShootingNewton.setProperty(100);
        RungeKuttaStages      = 2;
        StoreDecompositions   = false;
        TransientSolver       = false;
    end
    
    properties (SetAccess = protected)
        LinearSolves
        ShootingIterations
    end
    
    methods
        %% Constructor
        function this = ShootingNewton(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Setters
        function this = set.ShootingTolerance(this, tol)
            this.ShootingTolerance = this.setProperty(tol);
        end
        
        function this = set.NewtonTolerance(this, tol)
            this.NewtonTolerance = this.setProperty(tol);
        end
        
        function this = set.GMRESTolerance(this, tol)
            this.GMRESTolerance = this.setProperty(tol);
        end  
        
        function this = set.MaxShootingIterations(this, maxItt)
            this.MaxShootingIterations = this.setProperty(maxItt);
        end
        
        function this = set.MaxNewtonIterations(this, maxItt)
            this.MaxNewtonIterations = this.setProperty(maxItt);
        end
        
        function this = set.MaxGMRESIterations(this, maxItt)
            this.MaxGMRESIterations = this.setProperty(maxItt);
        end
        
        %% Solve
        function solution = solve(this, model, x0)
            %% Setup Matrices
            matrixFactory = DynamicMatrixFactory(copy(model));
            this.Matrices = matrixFactory;
            
            %% Get Algorithm Parameters
          	maxShootingItt = this.MaxShootingIterations.Value;
            maxNewtonItt   = this.MaxNewtonIterations.Value;
            maxGMRESItt    = this.MaxGMRESIterations.Value;
            shootingTol    = this.ShootingTolerance.Value;
            newtonTol      = this.NewtonTolerance.Value;
            gmresTol       = this.GMRESTolerance.Value;
            nStages        = this.RungeKuttaStages;
            [A,b,c]        = this.getButcherTable(nStages);
            storeLDL       = this.StoreDecompositions;
            verbose        = this.Verbose;
            
            %% Initialize
            times      = model.getTimePoints(this.TimePoints.Value);
            this.Times = times;
            nUnknowns  = length(matrixFactory.f(times(1)));   
            nTimes     = numel(times);
            x          = cell(1, nTimes);
            v          = cell(nStages, nTimes);
            d          = cell(1, nStages);
            x_t        = cell(1, nTimes);
            
            if nargin < 3 || isempty(x0)
                x{1} = zeros(nUnknowns, 1);
            else
                nAssemblies = numel(x0);
                x{1} = cell(nAssemblies, 1);
                for i = 1:nAssemblies
                    x{1}{i} = matrixFactory.PostProcessing(i).Full2Reduced * x0{i};
                end
                x{1} = cell2mat(x{1});
            end
            
            v{end,1} = zeros(nUnknowns, 1);
            v{end,1} = zeros(nUnknowns, 1);
            d{end}   = zeros(nUnknowns, 1);
            x_t{1}   = zeros(nUnknowns, 1);
            
            if storeLDL
                L = cell(1, nTimes*nStages);
                D = cell(1, nTimes*nStages);
                S = cell(1, nTimes*nStages);
                P = cell(1, nTimes*nStages);
                R = cell(1, nTimes*nStages);
            end
            
            iLinear     = 0;
            iShooting   = 1;
            shootingErr = 1;
            
            %% Begin Simulation Time
            tic
            while shootingErr > shootingTol && iShooting <= maxShootingItt
                for i = 2:nTimes
                    if verbose
                        sprintf('Shooting Iteration %d, Time %d', iShooting, i)
                    end
                        
                    %% For each time
                    h    = times(i) - times(i - 1);
                    x{i} = x{i - 1};
                    for s = 1:nStages
                        %% Calculate stage times
                        t  = times(i - 1) + h * c(s);
                        ha = h * A(s, s);
                        
                        %% Get initial guess
                        if iShooting == 1
                            v{s,i} = x{i - 1} + ha * x_t{i-1};
                        elseif this.TransientSolver
                            v{s,i} = 0 * x{i - 1};
                        end
                        
                        %% Create part of the stage derivative estimate
                        dPart = x{i - 1};
                        for j = 1:(s - 1);
                            dPart = dPart + A(s,j) * d{j};
                        end
                        
                        %% Create constant matrices
                        f  = matrixFactory.f(t, ha);
                        K  = matrixFactory.K(t, ha);
                        C  = matrixFactory.C(t, ha);
                        
                        %% Perform newton - raphson itteration
                        iNewton   = 1;
                        newtonErr = 1;
                        while newtonErr > newtonTol && iNewton < maxNewtonItt
                            [G,g]     = matrixFactory.G(t, v{s,i});
                            
                            J         = K + C + G;
                            
                            d{s}      = v{s, i} - dPart;
                            r         = C * d{s} + K * v{s,i} + g - f;
                            dv        = J \ r;
                            v{s,i}    = v{s,i} - dv;
                            
                            iNewton   = iNewton + 1;
                            newtonErr = norm(r) / norm(g - f);
                        end
                        
                        iLinear = iLinear + iNewton;
                        if storeLDL
                            [L{i}, D{i}, P{i}, S{i}] = ldl(J);
                            R{i} = -(K + G) / h / A(s, s);
                            iLinear = iLinear + 1;
                        end
                        
                        d{s} = d{s} / A(s,s);
                    end
                    
                    x{i}   = v{end,i};
                    x_t{i} = d{end} / h;
                end
                
                x_t{1} = x_t{end};
                for s = 1:nStages
                    v{s, 1} = v{s, end};
                end
                
                r  = x{end} - x{1};
                sf = sum(abs(J), 2);
                shootingErr = norm(r .* sf) / norm(x{end} .* sf);
                
                this.ConvergenceHistory(end+1) = shootingErr;
                
                if verbose
                    ShootingResidual = shootingErr
                end

                if shootingErr > shootingTol && maxGMRESItt > 0 && ~this.TransientSolver
                    if storeLDL
                        f = @(y)(ShootingNewton.NSTFJMVP(v,times,matrixFactory,A,b,c,y,L,D,P,S,R));
                        if verbose
                            dx = gmres(f, r, maxGMRESItt, gmresTol, 1, [], [], r);
                        else
                            [dx,~,~,~] = gmres(f, r, maxGMRESItt, gmresTol, 1, [], [], r);
                        end
                    else
                        f = @(y)(ShootingNewton.NSTFJMVP(v,times,matrixFactory,A,b,c,y));
                        [dx,~,~,itter] = gmres(f, r, maxGMRESItt, gmresTol, 1, [], [], r);
                        iLinear = iLinear + itter(1) * itter(2) * nStages * nTimes;
                    end
                    x{1} = x{1} - dx;
                else
                    x{1} = x{end};
                end
                iShooting = iShooting + 1;
            end
            this.LinearSolves       = iLinear;
            this.ShootingIterations = iShooting;
            
            %% End Simulation Time
            this.SimulationTime = toc;
            
            if verbose
                SimulationTime = this.SimulationTime
            end
            
            %% Save the solution
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x,x_t);
            solution           = Solution(this);
        end
    end
    
    methods (Static)
        function mvp = NSTFJMVP(v,times,matrixFactory,A,b,c,r,L,D,P,S,R)
            %% Nonlinear State Transition Function Jacobian Matrix-Vector Product
            if nargin <= 7
                matricesAreStored = false;
            else
                matricesAreStored = true;
            end
            
            nTimes  = length(times);
            nStages = length(b);
            w       = cell(nStages, nTimes);
            y       = cell(1, nTimes);
            y{1}    = r;
            for i = 2:nTimes
                h = times(i) - times(i - 1);
                y{i} = y{i - 1};
                for s = 1:nStages         
                    rhs = y{i - 1};
                    for j = 1:(s - 1);
                        rhs = rhs + h * A(s, j) * w{j,i};
                    end
                    
                    if matricesAreStored
                        w{s,i} = S{i}*P{i}*(L{i}'\(D{i}\(L{i}\(P{i}'*S{i}*R{i}*rhs))));
                    else
                        t      = times(i - 1) + h * c(s);
                        ha     = h * A(s, s);
                        K      = matrixFactory.K(t, ha);
                        C      = matrixFactory.C(t, ha);
                        G      = matrixFactory.G(t, v{s, i});
                        J      = K + C + G;
                        rhs    = -((K + G) * rhs) / h / A(s, s);
                        w{s,i} = J \ rhs;
                    end
                    
                    y{i}   = y{i} + h * b(s) * w{s,i};
                end
            end
            mvp = y{end} - y{1};
        end
        
        function [A,b,c] = getButcherTable(nStages)
            switch nStages
                case 1
                    A = 1;
                case 2
                    A = [6 - 3*sqrt(2), 0;
                         3*sqrt(2),     6 - 3*sqrt(2)] / 6;
                case 3
                    l = 0.435866521508459;
                    A = [l,                  0,               0;
                         (1-l)/2,            l,               0;
                         (-6*l^2+16*l-1)/4, (6*l^2-20*l+5)/4, l];
                case 4
                    A = [6 - 3*sqrt(2),   0,               0,           0;
                         3 * sqrt(2),     6 - 3 * sqrt(2), 0,           0;
                         30 - 18*sqrt(2), 24 * sqrt(2)-36, 6-3*sqrt(2), 0;
                         2 * sqrt(2) + 1, sqrt(2) - 2,     1,           6-3*sqrt(2)] / 6;
            end
            b = A(end, :);
            c = sum(A, 2);
        end
      
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = ShootingNewton(varargin{:});
            end
        end
    end
end