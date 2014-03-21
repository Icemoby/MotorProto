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
        GMRESTolerance        = ShootingNewton.setProperty(1e-6);
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
          	maxShootingIter = this.MaxShootingIterations.Value;
            maxNewtonIter   = this.MaxNewtonIterations.Value;
            maxGMRESIter    = this.MaxGMRESIterations.Value;
            shootingTol    = this.ShootingTolerance.Value;
            newtonTol      = this.NewtonTolerance.Value;
            gmresTol       = this.GMRESTolerance.Value;
            nStages        = this.RungeKuttaStages;
            [A,b,c]        = this.getButcherTable(nStages);
            
            alpha = inv(A);
            gamma = sum(alpha,2);
            
            %% Initialize
            times      = model.getTimePoints(this.TimePoints.Value);
            this.Times = times;
            nUnknowns  = length(matrixFactory.f(times(1)));   
            nTimes     = numel(times) - 1;
            
            x          = cell(1, nTimes+1);
            x_t        = cell(1, nTimes+1);
            y          = cell(nStages, nTimes+1);
            
            if nargin < 3 || isempty(x0)
                x{1} = zeros(nUnknowns, 1);
                y{end,1} = zeros(nUnknowns, 1);
                d_ik = 0 * y{end,1};
            else
                nAssemblies = numel(x0);
                x{1} = cell(nAssemblies, 1);
                for i = 1:nAssemblies
                    x{1}{i} = matrixFactory.PostProcessing(i).Full2Reduced * x0{i};
                end
                x{1} = cell2mat(x{1});
                y{end,1} = x{1};
                d_ik = 0 * y{end,1};%d_ik = x_t{1}
            end
            
            if this.StoreDecompositions
                Ca = cell(nStages, nStages, nTimes);
                Cg = cell(nStages, nTimes);
                L = cell(nStages, nTimes);
                D = cell(nStages, nTimes);
                S = cell(nStages, nTimes);
                P = cell(nStages, nTimes);
            end
            
            iLinear     = 0;
            iShooting   = 0;
            shootingRes = 1;
            
            %%Begin Simulation Timing
            if this.Verbose
                display(sprintf('Shooting-Newton %d/%d\n',nStages,nTimes));
            end
            
            tic
            while shootingRes > shootingTol && iShooting <= maxShootingIter
                iShooting = iShooting + 1;
                for k = 1:nTimes
                    hk  = times(k+1) - times(k);
                    for i = 1:nStages
                        %% Calculate stage times
                        t_ik = times(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        %% Get initial guess
                        if i == 1
                            y{1, k+1} = y{nStages, k} + hk * c(i) * d_ik;
                        else
                            y{i, k+1} = y{i-1, k+1} + hk * (c(i) - c(i-1)) * d_ik;
                        end
                        
                        for j = 1:i
                            Ca{i,j,k} = matrixFactory.C(t_ik, hk / alpha(i,j), h_ik);
                        end
                        Cg{i,k} = matrixFactory.C(t_ik, hk / gamma(i), h_ik);
                        
                        f_ik = matrixFactory.f(t_ik, h_ik);
                        K_ik = matrixFactory.K(t_ik, h_ik);
                        
                        %% Perform newton - raphson itteration
                        iNewton   = 0;
                        newtonErr = 1;
                        while newtonErr > newtonTol && iNewton < maxNewtonIter
                            iNewton = iNewton + 1;
                            [G_ik, g_ik] = matrixFactory.G(t_ik, y{i,k+1});
                            J_ik = K_ik + Ca{i,i,k} + G_ik;

                            r_ik = K_ik * y{i,k+1} + g_ik - f_ik;
                            for j = 1:i
                                r_ik = r_ik + Ca{i,j,k} * y{j,k+1};
                            end
                            r_ik = r_ik - Cg{i,k} * y{nStages, k};
                            
                            if this.StoreDecompositions && ~this.TransientSolver
                                [L_ik, D_ik, P_ik, S_ik] = ldl(J_ik);
                                r_ik = S_ik  * r_ik;
                                r_ik = P_ik' * r_ik;
                                r_ik = L_ik  \ r_ik;
                                r_ik = D_ik  \ r_ik;
                                r_ik = L_ik' \ r_ik;
                                r_ik = P_ik  * r_ik;
                                r_ik = S_ik  * r_ik;
                                y{i,k+1} = y{i,k+1} - r_ik;
                            else
                                r_ik = J_ik \ r_ik;
                                y{i,k+1} = y{i,k+1} - r_ik;
                            end
                            
                            newtonErr = norm(r_ik) / norm(g_ik - f_ik);
                        end
                        
                        d_ik = (-gamma(i) / hk) * y{end,k};
                        for j = 1:i
                            d_ik = d_ik + (alpha(i,j) / hk) * y{j,k+1};
                        end
                        
                        iLinear = iLinear + iNewton;
                        if this.StoreDecompositions && ~this.TransientSolver
                            L{i,k} = L_ik;
                            D{i,k} = D_ik;
                            P{i,k} = P_ik;
                            S{i,k} = S_ik;
                        end
                    end
                end
                
                r = y{end,end}-y{end,1};
                shootingRes = norm(r) / norm(y{end,end});
                
                this.ConvergenceHistory(end+1) = shootingRes;
                
                if this.Verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', iShooting, shootingRes));
                end

                if shootingRes > shootingTol && ~this.TransientSolver
                    if this.StoreDecompositions
                        f = @ShootingNewton.LDLStoredMVP;
                        if this.Verbose
                            dy = gmres(f, r, maxGMRESIter, gmresTol, 1, [], [], r, Ca, Cg, L, D, P, S);
                            display(sprintf(' '));
                        else
                            [dy,~,~,~] = gmres(f, r, maxGMRESIter, gmresTol, 1, [], [], r, Ca, Cg, L, D, P, S);
                        end
                    else
                        f = @(z)(ShootingNewton.UnstoredMVP(z,y,times,alpha,gamma,c,matrixFactory));
                        [dy,~,~,iter] = gmres(f, r, maxGMRESIter, gmresTol, 1, [], [], r);
                        iLinear = iLinear + iter(1) * iter(2) * nStages * nTimes;
                    end
                    y{end,1} = y{end,1} - dy;
                elseif this.TransientSolver
                    y{end,1} = y{end,end};
                    display(sprintf(' '));
                end
            end
            this.LinearSolves       = iLinear;
            this.ShootingIterations = iShooting;            
            
            for k = 1:nTimes
                hk = times(k+1) - times(k);
                
                x{k} = y{end,k};
                
                x_t{k+1} = (- gamma(end) / hk) * y{end,k};
                for i = 1:nStages
                    x_t{k+1} = x_t{k+1} + (alpha(end,i) / hk) * y{i,k+1};
                end
            end
            x{end} = x{1};
            x_t{1} = x_t{end};
            
            %% End Simulation Time
            this.SimulationTime = toc;
            
            if this.Verbose
                display(sprintf('Simulation Time = %0.3g seconds\n', this.SimulationTime));
            end
            
            %% Save the solution
%             if this.StoreDecompositions
%                 clearvars Ca Cg L D P S
%             end
            
            [x, x_t] = matrixFactory.doPostProcessing(x, x_t);
            this.X   = x;
            this.X_t = x_t;
            solution = Solution(this);
        end
    end
    
    methods (Static)
        function mvp = LDLStoredMVP(z_s0,Ca,Cg,L,D,P,S)
            [nStages, nTimes] = size(Cg);
            
            z_sk = z_s0;
            z = cell(nStages, 1);
            for k = 1:nTimes
                for i = 1:nStages         
                    rhs = Cg{i,k} * z_sk;
                    for j = 1:(i - 1);
                        rhs = rhs - Ca{i,j,k} * z{j};
                    end
                    z{i} = S{i,k}*(P{i,k}*(L{i,k}'\(D{i,k}\(L{i,k}\(P{i,k}'*(S{i,k}*rhs))))));
                end
                z_sk = z{end};
            end
            mvp = z_sk - z_s0;
        end
        
        function mvp = UnstoredMVP(z_s0,y,times,alpha,gamma,c,matrixFactory)
            [nStages, nTimes] = size(y);
            nTimes = nTimes - 1;
            
            Ca = cell(nStages, 1);
            z = cell(nStages, nTimes+1);
            z{end,1} = z_s0;
            for k = 1:nTimes
                hk = times(k+1) - times(k);
                for i = 1:nStages
                    t_ik = times(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    
                    for j = 1:i
                        Ca{j} = matrixFactory.C(t_ik, hk / alpha(i,j), h_ik);
                    end
                    Cg = matrixFactory.C(t_ik, hk / gamma(i), h_ik);
                    
                    r_ik = Cg * z{end, k};
                    for j = 1:(i-1)
                        r_ik = r_ik - Ca{j} * z{j,k+1};
                    end

                  	K_ik = matrixFactory.K(t_ik, h_ik);
                    G_ik = matrixFactory.G(t_ik, y{i,k+1});
                    J_ik = K_ik + Ca{i} + G_ik;
                    
                    z{i,k+1} = J_ik \ r_ik;
                end
            end
            mvp = z{end,end} - z{end,1};
        end
        
        function [A,b,c] = getButcherTable(nStages)
            switch nStages
                case 1
                    A = 1;
                case 2
                    A = [4 0;
                         9 3]/12;
                case 4
                  	A = [2  0 0 0;
                         3  3 0 0;
                         2 -2 4 0;
                         0  0 9 3]/12;
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