classdef ExplicitErrorCorrection < Solver
    properties
        ShootingTolerance      = ShootingNewton.setProperty(sqrt(eps));
        NewtonTolerance        = ShootingNewton.setProperty(sqrt(eps));
        GMRESTolerance         = ShootingNewton.setProperty(sqrt(eps));
        MaxShootingIterations = ShootingNewton.setProperty(1000);
        MaxNewtonIterations   = ShootingNewton.setProperty(100);
        MaxGMRESIterations    = ShootingNewton.setProperty(1000);
        RungeKuttaStages       = 3;
        StoreDecompositions    = false;
    end
    
    properties (SetAccess = protected)
        SimulationTime
        LinearSolves
        ShootingIterations
    end
    
    methods
        function this = ExplicitErrorCorrection(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        function this = set.ShootingTolerance(this,tol)
            this.Tolerance = this.setProperty(tol);
        end
        
        function this = set.NewtonTolerance(this,tol)
            this.NewtonTolerance = this.setProperty(tol);
        end
        
        function this = set.GMRESTolerance(this,tol)
            this.GMRESTolerance = this.setProperty(tol);
        end  
        
        function this = set.MaxShootingIterations(this,maxItt)
            this.MaxShootingIterations = this.setProperty(maxItt);
        end
        
        function this = set.MaxNewtonIterations(this,maxItt)
            this.MaxNewtonIterations = this.setProperty(maxItt);
        end
        
        function this = set.MaxGMRESIterations(this,maxItt)
            this.MaxGMRESIterations = this.setProperty(maxItt);
        end
        
        function solution = solve(this, mesh, x0)
            %% Setup Matrices
            matrixFactory = DynamicMatrixFactory(mesh);
            this.Matrices = matrixFactory;
            
            %% Get Algorithm Parameters
          	maxShootingItt = this.MaxShootingIterations.Value;
            maxNewtonItt   = this.MaxNewtonIterations.Value;
            shootingTol    = this.ShootingTolerance.Value;
            newtonTol      = this.NewtonTolerance.Value;
            
            %% Initialize
            times      = matrixFactory.getTimePoints(this.TimePoints.Value);
            this.Times = times;
            nUnknowns  = length(matrixFactory.f(times(1)));   
            nTimes     = numel(times);
            x          = cell(1,nTimes);
            x_t        = cell(1,nTimes);
            
            if nargin < 3 || isempty(x0)
                x{1} = zeros(nUnknowns,1);
            else
                nAssemblies = numel(x0);
                x{1} = cell(nAssemblies,1);
                for i = 1:nAssemblies
                    x{1}{i} = matrixFactory.PostProcessing(i).Full2Reduced * x0{i};
                end
                x{1} = cell2mat(x{1});
            end
            
            x_t{1}   = zeros(nUnknowns,1);
            
            iLinear     = 0;
            iShooting   = 1;
            shootingErr = 1;
            tic
            while shootingErr > shootingTol && iShooting <= maxShootingItt
                E      = cell(3,3);
                E{1,1} = sparse(nUnknowns,nUnknowns);
                E{1,2} = E{1,1};
                E{2,1} = E{1,1};
                E{2,2} = E{1,1};

                for i = 2:nTimes
                    %% For each time
                    h    = times(i) - times(i - 1);
                    
                    %% Calculate stage times
                    t  = times(i);
                    
                    %% Get initial guess
                    if iShooting == 1
                        x{i} = x{i - 1} + h * x_t{i - 1};
                    end
                    
                    %% Create constant matrices
                    f  = matrixFactory.f(t, h);
                    K  = matrixFactory.K(t, h);
                    C  = matrixFactory.C(t, h);
                        
                    %% Perform newton - raphson itteration
                    iNewton   = 1;
                    newtonErr = 1;
                    while newtonErr > newtonTol && iNewton < maxNewtonItt
                        [G,g]   = matrixFactory.G(t, x{i});
                        
                        J       = K + C + G;
                        
                        x_t{i}  = x{i} - x{i-1};
                        r       = C * x_t{i} + K * x{i} + g - f;
                        dx      = J \ r;
                        x{i}    = x{i} - dx;
                        
                        iNewton   = iNewton + 1;
                        newtonErr = norm(r) / norm(g - f);
                    end
                    E{1,1} = E{1,1} +           (K+G);
                    E{1,2} = E{1,2} + (i-1)   * (K+G) * h;
                    E{2,1} = E{2,1} + (i-1)   * (K+G) * h;
                    E{2,2} = E{2,2} + (i-1)^2 * (K+G) * h^2;
                    
                    iLinear = iLinear + iNewton;
                    x_t{i}  = x_t{i} / h;
                end
                f      = matrixFactory.f(times(2),h);
                K      = matrixFactory.K(times(2),h);
                [G,g]  = matrixFactory.G(times(2),x{2});
                C      = matrixFactory.C(times(2),h);
                E{2,2} = E{2,2} + h^2 * (nTimes - 2) * (nTimes - 1) * C / 2;
                
                r      = K * x{2} + C * (x{2} - x{end}) + g - f;
                shootingErr = norm(r) / norm(g-f)
                
                r      = cell(2,1);
                r{1}   = C * (x{end} - x{1});
                r{2}   = h * r{1};
                
                r      = cell2mat(r(1));
                E      = cell2mat(E(1,1));
                dx     = -E \ r;
                dx     = mat2cell(dx, repmat(numel(dx),1,1)/1, 1);
                rold   = r;
                for i = 2:nTimes
                    x{i} = x{i} + dx{1}*0;% + h * (i - 1) * dx{2};
                end
                x{1}      = x{end};
                x_t{1}    = (x{end} - x{end-1}) / h;
                iShooting = iShooting + 1;
            end
            
            this.SimulationTime     = toc;
            this.LinearSolves       = iLinear;
            this.ShootingIterations = iShooting;
            
            %% Save the solution
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x,x_t);
            solution           = Solution(this);
        end
    end
    
    methods (Static)
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = ExplicitErrorCorrection(varargin{:});
            end
        end
    end
end