classdef HarmonicBalance < Solver    
    %HarmonicBalance.m Frequency domain periodic steady-state solver
    %
    % HarmonicBalance properties:
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
        MaxNewtonIterations = 100;
        NewtonTolerance     = sqrt(eps);
    end
    
    methods
        function this = HarmonicBalance(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        function solution = solve(this, model, x0)
            %% Configure matrices
            matrixFactory = HarmonicMatrixFactory(copy(model));
            this.Matrices = matrixFactory;
            
            %% Configure algorithm
          	maxNIt = this.MaxNewtonIterations;
            resTol = this.NewtonTolerance;
            
            %% Get harmonics
            [t,h]      = matrixFactory.getTimePoints(this.TimePoints);
            this.Times = t;
            t(end)     = [];
            
            %% Begin Simulation Time
            tic
            
            %% Create constant matrices
            f = matrixFactory.f(t,h);
            H = matrixFactory.K(t,h) + matrixFactory.C(t,h);
            
            %% Get initial guess
            if nargin < 3
                x = zeros(size(f));
            else
                x = x0;
            end
            
            %% Perform newton - raphson itteration
            iIter   = 1;
            relRes  = 1;
            verbose = this.Verbose;
            while iIter <= maxNIt && relRes > resTol
                [G,g] = matrixFactory.G(x, t, h);

                r  = H * x + g - f;
                J  = H + G;
            	dx = J \ r;
              	x  = x - dx;

                iIter   = iIter + 1;
                relRes = norm(r) / norm(g - f);
                
                if verbose
                    HarmonicResidual = relRes
                end
                
                this.ConvergenceHistory(end+1) = relRes;
            end
            
            %% End Simulation Time
            this.SimulationTime = toc;
            this.SimulationTime
            
            %% Post Processing
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x, t, h);
            solution           = Solution(this);
        end
    end
    
    methods (Static)
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = HarmonicBalance(varargin{:});
            end
        end
    end
end