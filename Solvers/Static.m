classdef Static < Solver
    %Static.m The standard solver for magnetostatic problems.
    %   Static is a solver for periodic magnetostatic problems. For a given
    %   model, a series of magnetostatic simulations are performed at a
    %   number of time points over the problem period.
    %
    % Static properties:
    %   NewtonTolerance     - Sets the relative residual tolerance
    %   MaxNewtonIterations - Sets the maximum number of iterations
    %
    % Static inherits properties and methods Solver.
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
        NewtonTolerance     = Static.setProperty(sqrt(eps));
        MaxNewtonIterations = Static.setProperty(100);
    end
    
    methods
        function this = Static(varargin)
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
        
        function solution = solve(this, model, x0)
            %See also Solver/solve

            %% Configure matrices
            matrixFactory = StaticMatrixFactory(model);
            this.Matrices = matrixFactory;
            
            %% Configure algorithm
          	maxNIt  = this.MaxNewtonIterations.Value;
            nrrTol  = this.NewtonTolerance.Value;
            verbose = this.Verbose;
            
            %% Get times points
            t = model.getTimePoints(this.TimePoints.Value);
            
            if this.TimePoints.Value == 1
                t = t(1);
            end

            this.Times = t;
            Nt         = numel(t);
            x          = cell(1, Nt);
            
            %% For each time point, 
            tic
            
            for i = 1:Nt
                %% Create constant matrices,
                f = matrixFactory.f(t(i));
                K = matrixFactory.K(t(i));

                %% Get initial guess,
                if i > 1
                    x{i} = x{i-1};
                else
                    if nargin == 3
                        x{1} = x0;
                    else
                        x{i} = zeros(size(f));
                    end
                end
                
                %% Perform Newton-Raphson iteration, 
                nIter  = 1;
                relRes = 1;
                
                if verbose
                    Time = i
                end
                
                while nIter < maxNIt && relRes > nrrTol
                    [G, g] = matrixFactory.G(t(i), x{i});
                    
                    r = K * x{i} + g - f;
                    J = K + G;
                    
                   	if i == 1 && nIter == 1
                        [~,~,p,S] = ldl(J, 'vector');
                        [~, q]    = sort(p);
                        spparms('autoamd', 1);
                    end
                    
                    %Explicit symmetrization ensures MA57 is called for ldl decomposition
                    %SYMAMD preording is precalculated since all matrices have
                    %the same structure. spparms('autoamd',0) turns of this
                    %stage in the sparse solver.
                    
                    J  = (J + J.');
                    r  = 2 * r(p);
                    J  = J(p,p);
                    dx = J \ r;
                    dx = dx(q);
                    
                    %Equivalent direct call to ldl. Tends to be slower than the explicit
                    %symmetrization above. May be faster if (J+J.') is
                    %expensive to compute.

%                     r       = r(p);
%                     J       = J(p,p);
%                     [L,D,P] = ldl(J);
%                     dx      = (P' \ (L' \ (D \ (L \ (P \ r)))));
%                     dx      = dx(q);
                    
                    x{i} = x{i} - dx;
                    
                    nIter  = nIter + 1;
                    relRes = norm(S * r) / norm(S * (g - f));
                    
                    this.ConvergenceHistory(end+1) = relRes;
                end
            end
            spparms('autoamd', 1);
            
            this.SimulationTime = toc;
            
            if verbose
                SimulationTime = this.SimulationTime
            end
            
            %% Save Solution
            this.X = x;
            fe     = 1 / t(end);
            
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x, fe);
            solution           = Solution(this);
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
                solverOut = Static(varargin{:});
            end
        end
    end
end