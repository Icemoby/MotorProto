classdef TPFEM < Solver
    %TPFEM.m Time-Periodic FEM steady-state solver
    %
    % TPFEM properties:
    % 	NewtonTolerance       - Tolerance of iteration occuring at each time point
    % 	GMRESTolerance        - Tolerance of the initial condition correction
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
        NewtonTolerance       = 1e-6;
        GMRESTolerance        = 1e-6;
        MaxNewtonIterations   = 100;
        MaxGMRESIterations    = 100;
        RungeKuttaStages      = 2;
        StoreDecompositions   = true;
    end
    
    properties (SetAccess = protected)
        LinearSolves
        Iterations
    end
    
    methods
        %% Constructor
        function this = TPFEM(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Solve
        function solution = solve(this, model, x0)
            %% Setup Matrices
            matrixFactory = DynamicMatrixFactory(copy(model));
            this.Matrices = matrixFactory;
            
            %% Get Algorithm Parameters
            maxNewtonIter = this.MaxNewtonIterations;
            maxGMRESIter  = this.MaxGMRESIterations;
            newtonTol     = this.NewtonTolerance;
            gmresTol      = this.GMRESTolerance;
            nStages       = this.RungeKuttaStages;
            [A,b,c]       = this.getButcherTable(nStages);
            
            alpha = inv(A);
            gamma = sum(alpha,2);
            
            %% Initialize
            times      = model.getTimePoints(this.TimePoints);
            this.Times = times;
            nUnknowns  = length(matrixFactory.f(times(1)));   
            nTimes     = numel(times) - 1;
            
            y = cell(nStages, nTimes);
            [y{:}] = deal(zeros(nUnknowns, 1));
            r = cell(nStages, nTimes);
            [r{:}] = deal(zeros(nUnknowns, 1));
            
            Ca = cell(nStages, nStages, nTimes);
            Cg = cell(nStages, nTimes);
            L = cell(nStages, nTimes);
            D = cell(nStages, nTimes);
            S = cell(nStages, nTimes);
            P = cell(nStages, nTimes);
            
            iter = 0;
            relres = 1;
            
            %% Begin Simulation Timing
            if this.Verbose
                display(sprintf('TPFEM %d/%d\n',nStages,nTimes));
            end
            
            tic
            while relres > newtonTol && iter < maxNewtonIter
                iter = iter + 1;
                normf = 0;
                normr = 0;
                for k = 1:nTimes
                    km1 = mod(k-2, nTimes) + 1;
                    hk = times(k+1) - times(k);
                    for i = 1:nStages
                        t_ik = times(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        for j = 1:i
                            Ca{i,j,k} = matrixFactory.C(t_ik, hk / alpha(i,j), h_ik);
                        end
                        Cg{i,k} = matrixFactory.C(t_ik, hk / gamma(i), h_ik);
                        
                        f_ik = matrixFactory.f(t_ik, h_ik);
                        K_ik = matrixFactory.K(t_ik, h_ik);

                        [G_ik, g_ik] = matrixFactory.G(t_ik, y{i,k});
                        J_ik = K_ik + Ca{i,i,k} + G_ik;

                        r{i,k} = K_ik * y{i,k} + g_ik - f_ik;
                        for j = 1:i
                            r{i,k} = r{i,k} + Ca{i,j,k} * y{j,k};
                        end
                        r{i,k} = r{i,k} - Cg{i,k} * y{nStages, km1};

                      	[L{i,k}, D{i,k}, P{i,k}, S{i,k}] = ldl(J_ik);

                        normf = normf + norm(f_ik)^2;
                        normr = normr + norm(r{i,k})^2;
                    end
                end
                relres = sqrt(normr / normf);

                if this.Verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', iter, relres));
                end
                
                if relres > newtonTol
                    A = @TPFEM.MVP1;
                    dy = reshape(r,[],1);
                    dy = cell2mat(dy);
                    if this.Verbose
                        dy = gmres(A, dy, maxGMRESIter, gmresTol, 1, [], [], dy, Ca, Cg, L, D, P, S);
                        display(sprintf(' '));
                    else
                        [dy,~,~,~] = gmres(A, dy, maxGMRESIter, gmresTol, 1, [], [], dy);
                    end
                    dy = this.PC1(dy,Ca,Cg,L,D,P,S);
                    dy = mat2cell(dy,nUnknowns*ones(1,nTimes*nStages));
                    dy = reshape(dy,nStages,nTimes);

                    for i = 1:nTimes
                        for j = 1:nStages
                            y{j,i} = y{j,i} - dy{j,i};
                        end
                    end
                end
            end
            
           	x = cell(nTimes,1);
            x_t = cell(nTimes,1);
            for i = 1:nTimes
                hk = times(i+1) - times(i);
                
                j = mod(i-2,nTimes)+1;
                x{i} = y{end,j};
                
                k = mod(j-2,nTimes)+1;
                x_t{i} = (- gamma(end) / hk) * y{end,k};
                for k = 1:nStages
                    x_t{i} = x_t{i} + (alpha(end,k) / hk) * y{k,j};
                end
            end
            x{nTimes+1}   = x{1};
            x_t{nTimes+1} = x_t{1};
            
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
        function y = MVP1(x,Ca,Cg,L,D,P,S)
            [nStages, nTimes] = size(Cg);
            nUnknowns = length(x) / (nStages * nTimes);
            
            x = mat2cell(x, nUnknowns * ones(1, nStages * nTimes));
            x = reshape(x, nStages, nTimes);
            y = x;
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:nStages
                x{j,k} = S{j,k}*(P{j,k}*(L{j,k}'\(D{j,k}\(L{j,k}\(P{j,k}'*(S{j,k}*x{j,k}))))));
                for i = (j+1):nStages
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:nTimes                
                for j = 1:nStages
                    x{j,k} = x{j,k} + Cg{j,k} * x{nStages,k-1};
                end
                
                for j = 1:nStages
                    x{j,k} = S{j,k}*(P{j,k}*(L{j,k}'\(D{j,k}\(L{j,k}\(P{j,k}'*(S{j,k}*x{j,k}))))));
                    for i = (j+1):nStages
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            %% Matrix contribution from upper-right hand corner block
            for j = 1:nStages
                y{j,1} = y{j,1} - Cg{j,1} * x{nStages,nTimes};
            end
            
            y = reshape(y, nStages * nTimes, 1);
            y = cell2mat(y);
        end
        
        function x = PC1(x,Ca,Cg,L,D,P,S)
            [nStages, nTimes] = size(Cg);
            nUnknowns = length(x) / (nStages * nTimes);
            
            x = mat2cell(x, nUnknowns * ones(1, nStages * nTimes));
            x = reshape(x, nStages, nTimes);
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:nStages
                x{j,k} = S{j,k}*(P{j,k}*(L{j,k}'\(D{j,k}\(L{j,k}\(P{j,k}'*(S{j,k}*x{j,k}))))));
                for i = (j+1):nStages
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:nTimes                
                for j = 1:nStages
                    x{j,k} = x{j,k} + Cg{j,k} * x{nStages,k-1};
                end
                
                for j = 1:nStages
                    x{j,k} = S{j,k}*(P{j,k}*(L{j,k}'\(D{j,k}\(L{j,k}\(P{j,k}'*(S{j,k}*x{j,k}))))));
                    for i = (j+1):nStages
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            x = reshape(x, nStages * nTimes, 1);
            x = cell2mat(x);
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
                solverOut = TPFEM(varargin{:});
            end
        end
    end
end