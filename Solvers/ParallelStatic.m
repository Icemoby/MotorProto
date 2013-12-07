classdef ParallelStatic < Solver
    properties
        NewtonTolerance     = Static.setProperty(sqrt(eps));
        MaxNewtonIterations = Static.setProperty(100);
    end
    
    methods
        function this = ParallelStatic(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        function this = set.NewtonTolerance(this,tol)
            this.NewtonTolerance = this.setProperty(tol);
        end  
        
        function this = set.MaxNewtonIterations(this, maxNIt)
            this.MaxNewtonIterations = this.setProperty(maxNIt);
        end
        
        function solution = solve(this, model, ~)
            %% Configure matrices
            this.Matrices = StaticMatrixFactory(model);
            this.Times    = this.Matrices.getTimePoints(this.TimePoints.Value);
           	if this.TimePoints.Value == 1
                this.Times = this.Times(1);
            end
            
            %% Configure algorithm
            maxNIt  = this.MaxNewtonIterations.Value;
            nrrTol  = this.NewtonTolerance.Value;
            verbose = this.Verbose;           
            
            %% Setup Parallel Job
            tic
            N = 4;
            newWorkerPool(N - 1);
            toc
            scatterData(this, 'this', N, 'N', maxNIt, 'maxNIt', nrrTol, 'nrrTol');
            toc
            
            commands    = cell(1,4);
            commands{1} = ['warning off ', '''', 'MotorProto:Verbose', ''''];
            commands{2} = 'tic;I = setup(this, N, WID);';
            commands{3} = 't = this.Times(I);';
            commands{4} = 'x_loc = newtonRaphson(this, t, maxNIt, nrrTol);toc';
            send('ToAll', commands);
            toc
            
            I      = setup(this, N, 0);
            t      = this.Times(I);
            Nt     = numel(t);
            L      = length(this.Matrices.f(0));
            x      = zeros(L, Nt);
            x(:,I) = newtonRaphson(this, t, maxNIt, nrrTol);
            toc
            
            for i = 2:N
                J      = setup(this, N, i - 1);
                x(:,J) = retrieve(i - 1,'x_loc');
            end
            this.SimulationTime = toc;
            
            %% Save Solution
            this.X = x;
            fe     = 1 / this.Times(end);
            
            [this.X, this.X_t]  = this.Matrices.doPostProcessing(this.X, fe);
            solution            = Solution(this);
        end
        
        function I = setup(this, N, WID)
            t  = this.Times;
            Nt = numel(t);
            I  = (ceil(Nt / N * WID)+1):ceil(Nt / N * (WID + 1));
        end
        
        function x = newtonRaphson(this, t, maxNIt, nrrTol)
            Nt = numel(t);
            L  = length(this.Matrices.f(0));
            x  = zeros(L, Nt);
            
           	%% For each time point,
            for i = 1:Nt
                j = i;
                %% Create constant matrices,
                f = this.Matrices.f(t(j));
                K = this.Matrices.K(t(j));
                
                %% Get initial guess,
                if i > 1
                    x(:,i) = x(:,i-1);
                end
                
                %% Perform Newton-Raphson iteration,
                nIter  = 1;
                relRes = 1;
                
                while nIter < maxNIt && relRes > nrrTol
                    [G, g] = this.Matrices.G(t(j), x(:,i));
                    
                    r = K * x(:,i) + g - f;
                    J = K + G;
                    
                    if i == 1 && nIter == 1
                        [~,~,p,S] = ldl(J, 'vector');
                        [~, q]    = sort(p);
                        spparms('autoamd', 1);
                    end
                    
                    %Explicit symmetrization ensures MA57 is called for ldl decomposition
                    %SYMAMD preording is precalculated since all matrices have the same
                    %structure. spparms('autoamd',0) turns of this stage in the sparse solver.
                    
                    J  = (J + J.');
                    r  = 2 * r(p);
                    J  = J(p,p);
                    dx = J \ r;
                    dx = dx(q);
                    
                    %Equivalent direct call to ldl. Tends to be slower than the explicit
                    %symmetrization above. May be faster if (J+J.') is expensive to compute.
                    
                    %                     r       = r(p);
                    %                     J       = J(p,p);
                    %                     [L,D,P] = ldl(J);
                    %                     dx      = (P' \ (L' \ (D \ (L \ (P \ r)))));
                    %                     dx      = dx(q);
                    
                    x(:,i) = x(:,i) - dx;
                    nIter  = nIter + 1;
                    relRes = norm(S * r) / norm(S * (g - f));
                end
            end
            spparms('autoamd', 1);
        end
        
        function this = updateFrequency(this, fNew)
            fOld       = 1 / this.Times.Value(end);
            this.Times = this.Times.Value * fOld / fNew;
            x_t        = this.X_t;
            N          = numel(x_t);
            for i = 1:N
                x_t{i} = x_t{i} * fNew / fOld;
            end
        end
    end
    
    methods (Static)        
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = ParallelStatic(varargin{:});
            end
        end
    end
end