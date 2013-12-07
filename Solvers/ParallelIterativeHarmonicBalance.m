classdef ParallelIterativeHarmonicBalance < Solver
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
        function this = ParallelIterativeHarmonicBalance(varargin)
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
            nTimes     = numel(times);
            
            %% Setup Parallel Job
            WID = 0;
            N   = 4;
            newWorkerPool(N - 1);
            scatterData(this, 'this', N, 'N');
            
            commands    = cell(1,1);
            commands{1} = '[I, t, M, x, x_t, g, f, r, K, G, L1, U1, P1, Q1, R1, L2, U2, P2, Q2, R2, C, Omega, w] = setup(this, N, WID);';
            send('ToAll', commands);
            
            %% Initialize
            [I, t, M, x, x_t, g, f, r, K, G, L1, U1, P1, Q1, R1, L2, U2, P2, Q2, R2, C, Omega, w] = setup(this, N, WID);
            
            iShooting   = 1;
            shootingErr = 1;
            
            %% Begin Simulation Time
            tic
            commands{1} = '[r, K, G, L1, U1, P1, Q1, R1, D] = decomposeTimeDomain(this, t, I, x, x_t, f, g, r, K, G, C, L1, U1, P1, Q1, R1);';
            commands{2} = '[L2, U2, P2, Q2, R2]             = decomposeFreqDomain(this, w, D, C, L2, U2, P2, Q2, R2);';
            while shootingErr > shootingTol && iShooting <= maxShootingItt
                tic
                send('ToAll', commands(1));
                toc
                
                [r, K, G, L1, U1, P1, Q1, R1, D] = decomposeTimeDomain(this, t, I, x, x_t, f, g, g, K, G, C, L1, U1, P1, Q1, R1);
                
                tic
                D                                = retrieveDiagonal(this, D, N, nTimes);
                toc
                
                tic
                scatterData(D, 'D');
                toc
                
                tic
                send('ToAll', commands(2));
                toc
                
                [L2, U2, P2, Q2, R2]             = decomposeFreqDomain(this, w, D, C, L2, U2, P2, Q2, R2);
                
                tic
                r                                = retrieveVariable(this, r, 'r', N);
                toc

                r                                = reshape(r,[],1);
                
                A                                = @(x)(AFun(this, x, nTimes, I, N, L1, U1, P1, Q1, R1, Omega, K, G, C));
                M                                = @(x)(MFun(this, x, nTimes, I, N, L1, U1, P1, Q1, R1, Omega, K, G, C));
                dx                               = gmres(A, r, [], [], 100);
                dx                               = M(dx);
                
                dx                               = reshape(dx,[],nTimes - 1);
                x                                = x - dx;
                x_f                              = fft(x, [], 2) / (nTimes - 1);
                x_t                              = x_f * Omega;
                x_t                              = ifft(x_t, [], 2) * (nTimes - 1);
                
                tic
                scatterData(x, 'x', x_f, 'x_f', x_t, 'x_t');
                toc
                
                shootingErr = max(max(abs(dx),[],2) ./ max(abs(x),[],2))
            end

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
        
        function I = getLocalIndices(this, N, WID)
            t  = this.Times;
            Nt = numel(t);
            T  = t(end);
            h  = Nt / 2 - 1;
            Nt = Nt - 1;
            I  = (ceil(Nt / N * WID)+1):ceil(Nt / N * (WID + 1));
            t  = t(I);
            M  = numel(t);
        end
        
        function [I, t, M, x, x_t, g, f, r, K, G, L1, U1, P1, Q1, R1, L2, U2, P2, Q2, R2, C, Omega, w] = setup(this, N, WID)
            t  = this.Times;
            Nt = numel(t);
            T  = t(end);
            h  = Nt / 2 - 1;
            Nt = Nt - 1;
            I  = (ceil(Nt / N * WID)+1):ceil(Nt / N * (WID + 1));
            t  = t(I);
            M  = numel(t);
            
            nUnknowns = length(this.Matrices.f(t(1)));
            
            x   = zeros(nUnknowns, Nt);
            x_t = zeros(nUnknowns, Nt);
            K   = cell(1, M);
            G   = cell(1, M);
            L1  = cell(1, M);
            U1  = cell(1, M);
            P1  = cell(1, M);
            Q1  = cell(1, M);
            R1  = cell(1, M);
            L2  = cell(1, M);
            U2  = cell(1, M);
            P2  = cell(1, M);
            Q2  = cell(1, M);
            R2  = cell(1, M);
            
            g   = zeros(nUnknowns, M);
            f   = zeros(nUnknowns, M);
            r   = zeros(nUnknowns, M);
            
            C     = this.Matrices.C(0, 1);
            Omega = 2 * pi / T * 1i * sparse(1:(2*h+1), 1:(2*h+1), [0:h -h:1:-1]);
            w     = diag(Omega);
            w     = w(I);
        end
        
        function [r, K, G, L1, U1, P1, Q1, R1, D] = decomposeTimeDomain(this, t, I, x, x_t, f, g, r, K, G, C, L1, U1, P1, Q1, R1)
            nTimes = numel(t);
            for i = 1:nTimes
                f(:,i)        = this.Matrices.f(t(i), 1);
                K{i}          = this.Matrices.K(t(i), 1);
                [G{i},g(:,i)] = this.Matrices.G(t(i), real(x(:,I(i))));
                
                if i == 1
                    D = K{i} + G{i};
                else
                    D = D + K{i} + G{i};
                end
                
                [L1{i}, U1{i}, P1{i}, Q1{i}, R1{i}] = lu(K{i} + G{i});
                
                r(:,i) = K{i} * x(:,I(i)) + C * x_t(:,I(i)) + g(:,i) - f(:,i);
            end
        end
        
        function [L2, U2, P2, Q2, R2] = decomposeFreqDomain(this, w, D, C, L2, U2, P2, Q2, R2)
            nFreqs = numel(w);
            for i = 1:nFreqs
                [L2{i}, U2{i}, P2{i}, Q2{i}, R2{i}] = lu(D + C*w(i));
            end
        end
        
        function b = AFun(this, x, nTimes, I, N, L1, U1, P1, Q1, R1, Omega, K, G, C)
            x   = reshape(x, [], nTimes - 1);
            x_f = fft(x, [], 2) / (nTimes - 1);
            x_t = x_f * Omega;
            x_t = ifft(x_t, [], 2) * (nTimes - 1);
            
            tic
            scatterData(x, 'x', x_t, 'x_t');
            toc
            
            %% First PC
            tic
            send('ToAll',{'y = localPCApp(this, L1, U1, P1, Q1, R1, x(:,I));'});
            toc
            
            y      = zeros(length(x), nTimes - 1);
            y(:,I) = localPCApp(this, L1, U1, P1, Q1, R1, x(:,I));
            y      = retrieveVariable(this, y, 'y', N);
            y_t    = fft(y, [], 2) * Omega;
            y_t    = ifft(y_t, [], 2);
            
            tic
            scatterData(y, 'y', y_t, 'y_t');
            toc
            
            %% RHS Update
            tic
            send('ToAll',{'z = localMVP(this, I, K, G, C, y, y_t);'});
            toc
            
            z      = zeros(length(x), nTimes - 1);
            z(:,I) = localMVP(this, I, K, G, C, y, y_t);
            
            tic
            z      = retrieveVariable(this, z, 'z', N);
            toc
            
            %% Exit
            b = reshape(z, [], 1);
        end
        
        function b = MFun(this, x, nTimes, I, N, L1, U1, P1, Q1, R1, Omega, K, G, C)
            x   = reshape(x, [], nTimes - 1);
            x_f = fft(x, [], 2) / (nTimes - 1);
            x_t = x_f * Omega;
            x_t = ifft(x_t, [], 2) * (nTimes - 1);
            scatterData(x, 'x', x_t, 'x_t');
            
            %% First PC
            send('ToAll',{'y = localPCApp(this, L1, U1, P1, Q1, R1, x(:,I));'});
            y      = zeros(length(x), nTimes - 1);
            y(:,I) = localPCApp(this, L1, U1, P1, Q1, R1, x(:,I));
            y      = retrieveVariable(this, y, 'y', N);
            
            b      = reshape(y, [], 1);
            
%             y_t    = fft(y, [], 2) * Omega;
%             y_t    = ifft(y_t, [], 2);
%             
%             
%             scatterData(y, 'y', y_t, 'y_t');
%             
%             %% RHS Update
%             send('ToAll',{'z = localMVP(this, I, K, G, C, y, y_t);'});
%             
%             z      = zeros(length(x), nTimes - 1);
%             z(:,I) = localMVP(this, I, K, G, C, y, y_t);
%             z      = retrieveVariable(this, z, 'z', N);
%             
%             %% Exit
%             b = reshape(z, [], 1);
        end
        
        function b = localMVP(this, I, K, G, C, x, x_t)
            b = x(:,I) * 0;
            N = numel(I);
            for i = 1:N
                j = I(i);
                b(:,i) = K{i} * x(:,j) + G{i} * x(:,j) + C * x_t(:,j);
            end
        end
        
        function x = localPCApp(this, L, U, P, Q, R, b)
            x = b * 0;
            N = numel(L);
            for i = 1:N
                x(:,i) = (Q{i} * (U{i} \ (L{i} \ (P{i} * (R{i} \ b(:,i))))));
            end
        end
        
        function y = retrieveVariable(this, y, varString, N)
            for i = 2:N
                I      = getLocalIndices(this, N, i-1);
                y(:,I) = retrieve(i - 1,varString);
            end
        end
        
        function D = retrieveDiagonal(this, D, N, nTimes)
            for i = 2:N
                D = D + retrieve(i - 1, 'D');
            end
            D = D / (nTimes - 1);
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
                solverOut = ParallelIterativeHarmonicBalance(varargin{:});
            end
        end
    end
end