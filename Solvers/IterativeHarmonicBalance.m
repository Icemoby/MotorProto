classdef IterativeHarmonicBalance < Solver
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
        ShootingTolerance     = ShootingNewton.setProperty(sqrt(eps));
        NewtonTolerance       = ShootingNewton.setProperty(sqrt(eps));
        GMRESTolerance        = ShootingNewton.setProperty(sqrt(eps));
        MaxShootingIterations = ShootingNewton.setProperty(100);
        MaxNewtonIterations   = ShootingNewton.setProperty(100);
        MaxGMRESIterations    = ShootingNewton.setProperty(100);
        RungeKuttaStages      = 2;
        StoreDecompositions   = false;
    end
    
    properties (SetAccess = protected)
        LinearSolves
        ShootingIterations
    end
    
    methods
        %% Constructor
        function this = IterativeHarmonicBalance(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Setters
        function this = set.ShootingTolerance(this, tol)
            this.Tolerance = this.setProperty(tol);
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
            shootingTol    = this.ShootingTolerance.Value;
            verbose        = this.Verbose;
            
            %% Initialize
            times      = matrixFactory.getTimePoints(this.TimePoints.Value);
            times      = linspace(0,times(end),78);
%             times      = linspace(0,times(end), numel(times)+1);
            this.Times = times;
            nUnknowns  = length(matrixFactory.f(times(1)));   
            nTimes     = numel(times);
            x          = zeros(nUnknowns, nTimes-1);
            x_t        = zeros(nUnknowns, nTimes-1);
            K          = cell(1, nTimes-1);
            G          = cell(1, nTimes-1);
            L1         = cell(1, nTimes-1);
            U1         = cell(1, nTimes-1);
            P1         = cell(1, nTimes-1);
            Q1         = cell(1, nTimes-1);
            R1         = cell(1, nTimes-1);
            L2         = cell(1, nTimes-1);
            U2         = cell(1, nTimes-1);
            P2         = cell(1, nTimes-1);
            Q2         = cell(1, nTimes-1);
            R2         = cell(1, nTimes-1);
            
            s          = cell(100,1);
            g          = zeros(nUnknowns, nTimes-1);
            f          = zeros(nUnknowns, nTimes-1);
            r          = zeros(nUnknowns, nTimes-1);
            
            iShooting   = 1;
            shootingErr = 1;
            
            %% Begin Simulation Time
            tic
            C  = matrixFactory.C(0, 1);
            If = (any(C,1)) | (any(C,2).');
            It = ~If;
            h     = nTimes/2 - 1;
            Omega = 2*pi / times(end) * 1i * sparse(1:(2*h+1), 1:(2*h+1), [0:h -h:1:-1]);
            w     = diag(Omega);
            while shootingErr > shootingTol && iShooting <= maxShootingItt
                parfor i = 1:(nTimes-1)
                    t             = times(i);
                    f(:,i)        = matrixFactory.f(t, 1);
                    K{i}          = matrixFactory.K(t, 1);
                    [G{i},g(:,i)] = matrixFactory.G(t, real(x(:,i)));
                    
                    [L1{i},U1{i},P1{i},Q1{i},R1{i}] = lu(K{i}+G{i});
                end
                
                parfor i = 1:(nTimes-1)
                    r(:,i) = K{i}*x(:,i)+ C*x_t(:,i) + g(:,i) - f(:,i);
                end
                
                D = K{1}+G{1};
                for i = 2:(nTimes-1)%parfor turns sparse reduction variables into doubles
                    D = D + K{i}+G{i};
                end
                D = D / (nTimes-1);
                
                parfor i = 1:(nTimes-1)
                    [L2{i},U2{i},P2{i},Q2{i},R2{i}] = lu(D + C * w(i));
                end
                
                r = fft(r, [], 2) / (nTimes - 1);
                g = fft(g, [], 2) / (nTimes - 1);
                f = fft(f, [], 2) / (nTimes - 1);
                
                r = reshape(r,[],1);
                g = reshape(g,[],1);
                f = reshape(f,[],1);
                
                shootingErr = norm(r) / norm(g - f)
                
                A = @(x)(this.mvp(K, G, C, Omega, x));
%                 M = @(x)(this.pc(K, G, C, Omega, isLinear, x));
%                 M = @(x)(this.pc2(K, G, C, Omega, x));
                M = @(x)(this.pcfast(K,G,C,Omega,L1,U1,P1,Q1,R1,L2,U2,P2,Q2,R2,x));
%                 M = @(x)(this.pcsplit(K,G,C,Omega,It,If,x));
%                 AM = @(x)(A(M(x)));
                if iShooting == 1
                    dx  = gmres(A, r, [], [], 200, [], M);
                    s{1} = dx;
                else
                    dx                = gmres(A, r, [], [], 200, M, [], dx);
                    s{1}(:,iShooting) = s{1}(:,iShooting-1) + dx;
                end
                
                dx = reshape(dx, [], (nTimes-1));
                dx = ifft(dx, [], 2, 'symmetric') * (nTimes - 1);
                
                x  = x - dx;
                dx = fft(dx, [], 2) / (nTimes - 1);
                dx = reshape(dx, [], 1);
                
                x_t = fft(x, [], 2) * Omega;
                x_t = ifft(x_t, [], 2, 'symmetric');
                
                r = reshape(r,[],(nTimes - 1));
                g = reshape(g,[],(nTimes - 1));
                f = reshape(g,[],(nTimes - 1));
                iShooting = iShooting + 1
            end
            toc
            %% End Simulation Time
            this.SimulationTime = toc
            
            if verbose
                SimulationTime = this.SimulationTime
            end
            
            x   = mat2cell(x,   length(x),   ones(1, (nTimes - 1)));
            x_t = mat2cell(x_t, length(x_t), ones(1, (nTimes - 1)));
            
            x{end+1}   = x{1};
            x_t{end+1} = x_t{1};
            %% Save the solution
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x,x_t);
            solution           = Solution(this);
        end
    end
    
    methods (Static)
        function b = mvp(K, G, C, Omega, x)
            n   = length(K);
            
            x   = reshape(x,[],n);
            x_t = x * Omega;
            
            x   = ifft(x,   [], 2, 'symmetric') * n;
            x_t = ifft(x_t, [], 2, 'symmetric') * n;
            
            b   = 0 * x;
            parfor i = 1:n
                b(:,i) = K{i}*x(:,i)+G{i}*x(:,i)+C*x_t(:,i);
            end
            
            b = fft(b, [], 2) / n;
            b = reshape(b, [], 1);
        end
        
        function b = pc(K, G, C, Omega, I, x)
            n = length(K);
            
            x = reshape(x,[],n);
            
            if any(I)
                D = K{1}(I,I)+G{1}(I,I);
                for i = 2:n
                    D = D + K{i}(I,I)+G{i}(I,I);
                end
                D = D / n;

                b = 0 * x;
                for i = 1:n
                    b(I,i) = (D + C(I,I) * Omega(i,i)) \ x(I,i);
                end
            end
            
            I = ~I;
            if any(I)
                x = ifft(x, [],2) * n;
                for i = 1:n
                    b(I,i) = (K{i}(I,I)+G{i}(I,I)) \ x(I,i);
                end
                b(I,:) = fft(b(I,:), [], 2) / n ;
            end
            b = reshape(b,[],1);
        end
        
        function x = pc2(K, G, C, Omega, b)
            n = length(K);
            
            b = reshape(b,[],n);
            x = 0 * b;
            
            %% Second Iteration
            b = ifft(b, [], 2, 'symmetric') * n;
            x = ifft(x, [], 2, 'symmetric') * n;
            
            for i = 1:n
                x(:,i) = x(:,i) + (K{i}+G{i}) \ b(:,i);
            end
            x = fft(x, [], 2) / n;
            b = fft(b, [], 2) / n;
            
            %% RHS Update
            b = reshape(b,[],1);
            x = reshape(x,[],1);
            
            b = b - IterativeHarmonicBalance.mvp(K, G, C, Omega, x);
            
            b = reshape(b,[],n);
            x = reshape(x,[],n);
            
            %% First Iteration
            D = K{1}+G{1};
            for i = 2:n
                D = D + K{i}+G{i};
            end
            D = D / n;
            
            for i = 1:n
                x(:,i) = x(:,i) + (D + C * Omega(i,i)) \ b(:,i);
            end
            
            %% Exit
            x = reshape(x,[],1);
        end
      
        function x = pcfast(K, G, C, Omega, L1, U1, P1, Q1, R1, L2, U2, P2, Q2, R2, b)
            n = length(K);
            
            b = reshape(b,[],n);
            x = 0 * b;
            
            %% First Iteration
            b = ifft(b, [], 2, 'symmetric') * n;
            x = ifft(x, [], 2, 'symmetric') * n;
            
            for i = 1:n
                x(:,i) = x(:,i) + (Q1{i} * (U1{i} \ (L1{i} \ (P1{i} * (R1{i} \ b(:,i))))));
            end
            
            x = fft(x, [], 2) / n;
            b = fft(b, [], 2) / n;
            
            %% RHS Update
            b = reshape(b,[],1);
            x = reshape(x,[],1);
            
            b = b - IterativeHarmonicBalance.mvp(K, G, C, Omega, x);
            
            b = reshape(b,[],n);
            x = reshape(x,[],n);
            
            %% Second Iteration
            for i = 1:n
                x(:,i) = x(:,i) + (Q2{i} * (U2{i} \ (L2{i} \ (P2{i} * (R2{i} \ b(:,i))))));
            end
            
            %% Exit
            x = reshape(x,[],1);
        end
        
        function x = pcsplit(K, G, C, Omega, It, If, b)
            tic
            n = length(K);
            w = 1;
            r = b;
            b = reshape(b,[],n);
            x = 0 * b;
            
          	%% Frequency-Domain Part
            D = K{1}+G{1};
            for i = 2:n
                D = D + K{i}+G{i};
            end
            D = D / n;
            
            parfor i = 1:n
                dx      = (D(If,If) + C(If,If) * Omega(i,i)) \ b(If,i);
                x(If,i) = x(If,i) + w * dx;
            end
            
            %% Time-Domain Part
            b = ifft(b, [], 2, 'symmetric') * n;
            x = ifft(x, [], 2, 'symmetric') * n;
            
            parfor i = 1:n
%                 b(It,i) = b(It,i) - K{i}(It,If) * x(If,i) - G{i}(It,If) * x(If,i);
                dx      = (K{i}(It,It)+G{i}(It,It)) \ b(It,i);
                x(It,i) = x(It,i) + w * dx;
            end
            x = fft(x, [], 2) / n;
            
%             %% RHS Update
%             x = reshape(x,[],1);
%             
%             b = r - IterativeHarmonicBalance.mvp(K, G, C, Omega, x);
%             norm(b)/norm(r)
%             
%             b = reshape(b,[],n);
%             x = reshape(x,[],n);
% 
%             %% Time-Domain Part
%             b = ifft(b, [], 2, 'symmetric') * n;
%             x = ifft(x, [], 2, 'symmetric') * n;
%             
%             for i = 1:n
%                 dx      = (K{i}(It,It)+G{i}(It,It)) \ b(It,i);
%                 x(It,i) = x(It,i) + w * dx;
%                 b(If,i) = b(If,i) - K{i}(If,It) * dx - G{i}(If,It) * dx;
%             end
%             x = fft(x, [], 2) / n;
%             b = fft(b, [], 2) / n;
%             
%             %% Frequency-Domain Part
%             for i = 1:n
%                 dx      = (D(If,If) + C(If,If) * Omega(i,i)) \ b(If,i);
%                 x(If,i) = x(If,i) + w * dx;
%             end
%             
            %% Exit
            x = reshape(x,[],1);
            
            b = r - IterativeHarmonicBalance.mvp(K, G, C, Omega, x);
            norm(b)/norm(r)
            
            toc
        end
        
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = IterativeHarmonicBalance(varargin{:});
            end
        end
    end
end