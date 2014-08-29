classdef HarmonicBalance < Solver
    methods
        %% Constructor
        function this = HarmonicBalance(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Solve
        function solution = solve(this, model, ~)
            tic
            %% Setup Matrices
            getMatrix = DynamicMatrixFactory(model);
            this.Matrices = getMatrix;
            Nx = length(getMatrix.f(0));
            
            %% Initialize
            t = model.getTimePoints(this.TimePoints);
            T = t(end);
            t(end) = [];
            Nt = numel(t);
            Nh = Nt / 2;
            Y = zeros(Nx, Nt);
            
            %% Set Parameters
            errTol = 1e-2;
            errorEstimate = 1;
            
            gmresTol = 1e-2;
            gmresIter = 10;
            
            minIter = 0;
            newtonIter = 20;
            newtonTol = 0.1;
            finalTol = 1e-4;

            final = false;
            verbose = true;
            
            APC = @(x,K,G,C,M1,M2,h)(HarmonicBalance.MVP_PC0_BDTF(x,K,G,C,M1,M2,h));
            PC = @(x,K,G,C,M1,M2,h)(HarmonicBalance.PC0_BDTF(x,K,G,C,M1,M2,h));
                
            %% Allocate
            C = getMatrix.C(0,1,1);
            
            K = cell(1,Nt);
            f = zeros(Nx,Nt);
            
            for i = 1:Nt
                K{i} = getMatrix.K(t(i),1);
                f(:,i) = getMatrix.f(t(i),1);
            end
            
            while (errorEstimate > errTol || ~final) || (errorEstimate < 2*errTol && ~final)
                if (errorEstimate > 2 * errTol)
                    hx = [0:(Nt/2-1), 0, -(Nt/2-1):1:-1] * 2 * pi / T;
                    hy = [0:(2*Nt/2-1), 0, -(2*Nt/2-1):1:-1] * 2 * pi / T;

                    X = [Y(:,1:Nh), zeros(Nx, Nt-2*Nh),  Y(:,(Nh+1):end)];
                    Y = [Y(:,1:Nh), zeros(Nx, 2*(Nt-Nh)),Y(:,(Nh+1):end)];
                    y = ifft(Y,[],2,'symmetric') * (2*Nt);

                    Nh = Nt;
                    Nt = 2*Nt;
                    
                    if verbose
                        display(sprintf('Time Points = %d\n',Nh));
                    end
                    
                    t = linspace(0,T,Nt+1);
                    
                    G = cell(1,Nh);
                    g = zeros(Nx,Nt);
                    
                    K = [K, cell(1, Nt/2)];
                    f = [f, zeros(Nx, Nt/2)];
                    
                    for i = (Nt-1):-2:1
                        j = (i+1) / 2;
                        K{i} = K{j};
                        f(:,i) = f(:,j);
                    end
                    
                    for i = 2:2:Nt
                        K{i} = getMatrix.K(t(i),1);
                        f(:,i) = getMatrix.f(t(i),1);
                    end
                    
                    s = zeros(Nx,Nt);
                    
                    M1 = cell(5,Nh);
                    M2 = cell(5,Nh);
                    
                    h2 = norm(hx)^2 / Nh;
                else
                    if verbose
                        display(sprintf('Colocating to a tolerance of %f with %d time points.\n', finalTol, Nh));
                    end
                    
                    final = true;
                end
                
                iter = 0;
                res = 1;
                tol = 0;
                
                %% Start Iteration
                while ((res > tol) && (iter < newtonIter)) || (iter < minIter)
                    %% Calculate Residual
                    for i = 1:Nt
                        if mod(i,2) == 1
                            j = (i+1) / 2;
                            [G{j},g(:,i)] = getMatrix.G(t(i), y(:,i), 1);
                            s(:,i) = K{i} * y(:,i) + g(:,i) - f(:,i);
                        else
                            [~, g(:,i)] = getMatrix.G(t(i), y(:,i), 1);
                            s(:,i) = K{i} * y(:,i) + g(:,i) - f(:,i);
                        end
                    end
                    r = s(:,1:2:end);

                    R = fft(r,[],2) / Nh;
                    for i = 1:Nh
                        R(:,i) = R(:,i) + 1i*hx(i)*(C*X(:,i));
                    end
                    r = ifft(R,[],2,'symmetric') * Nh;

                    S = fft(s,[],2) / Nt;
                    Y = fft(y,[],2) / Nt;
                    for i = 1:Nt
                        S(:,i) = S(:,i) + 1i*hy(i)*(C*Y(:,i));
                    end
                    s = ifft(S,[],2,'symmetric') * Nt;
                    
                    
                    if ~final
                        res = 0;
                        for i = 1:Nh
                            res = res + norm(R(:,i))^2;
                        end
                        res = sqrt(res);

                        tol = 0;
                        for i = 1:Nt
                            tol = tol + norm(S(:,i))^2;
                        end
                        tol = sqrt(tol) * newtonTol;
                    else
                        res = 0;
                        for i = 1:Nh
                            res = max(res, norm(r(:,i)) / norm(f(:,2*i-1)-g(:,2*i-1)));
                        end
                        tol = finalTol;
                    end
                    
                    iter = iter + 1;
                    
                    if verbose
                        display(sprintf('Iteration %d, Relative Residual = %f\n', iter, res / tol));
                    end
                    
                    %% Solve Linearized Equation
                    if (res > tol) || (iter <= minIter)
                       	J0 = K{1} + G{1};
                        for i = 2:Nh
                            j = 2*i -1;
                            J0 = J0 + K{j} + G{i};
                        end
                        J0 = J0 / Nh;

                        for i = 1:Nh
                            j = 2*i - 1;
                            M1{1,i} = K{j}+G{i};
                            M2{1,i} = J0 + 1i*hx(i)*C;

                            [M1{3,i},M1{4,i},M1{2,i},M1{5,i},M1{1,i}] = lu(M1{1,i});
                            [M2{3,i},M2{4,i},M2{2,i},M2{5,i},M2{1,i}] = lu(M2{1,i});
                        end
                        
                        dX = gmres(APC, r(:), gmresIter, gmresTol, 1, [], [], r(:), K, G, C, M1, M2, hx);
                        dX = PC(dX,K,G,C,M1,M2,hx);
                        X  = X - dX;

                        Y(:,1:(Nh/2))       = X(:,1:(Nh/2));
                        Y(:,Nh/2+1)         = X(:,Nh/2+1) / 2;
                        Y(:,Nt-Nh/2+1)      = conj(X(:,Nh/2+1)) / 2;
                        Y(:,(Nt-Nh/2+2):Nt) = X(:,(Nh/2+2):Nh);
                        y = ifft(Y,[],2,'symmetric') * Nt;
                    end
                end
                
                resvec = zeros(1,Nt);
                for i = 1:Nt
                    resvec(i) = norm(s(:,i))/norm(f(:,i)-g(:,i));
                end
                errorEstimate = max(resvec);
            
               	if verbose
                    display(sprintf('Error Estimate = %f\n', errorEstimate));
                end
                
                figure(1);clf;bar(resvec)
                pause(0.1);
            end
            this.SimulationTime = toc;
            
           	if this.Verbose
                display(sprintf('Simulation Time = %0.3g seconds\n', this.SimulationTime));
            end
            
            %% Post Processing
            x = ifft(X, [], 2, 'symmetric') * Nh;
            x_t = X;
            for i = 1:Nh
                x_t(:,i) = 1i*hx(i)*x_t(:,i);
            end
            x_t = ifft(x_t, [], 2, 'symmetric') * Nh;
            
            x = [x, x(:,1)];
            x_t = [x_t, x_t(:,1)];
            
            x = mat2cell(x, Nx, ones(1,Nh+1));
            x_t = mat2cell(x_t, Nx, ones(1,Nh+1));
            
            [x, x_t] = getMatrix.doPostProcessing(x, x_t);
            this.Times = linspace(0,T,Nh+1);
            this.Y   = y;
            this.X   = x;
            this.X_t = x_t;
            
            solution = Solution(this);
        end
    end
    
    methods (Static)
    	function x = MVP_PC0_BDTF(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            z = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                %y(:,i) = M1{i} \ x(:,i); % Apply first preconditioner, y <- x1
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
            end
            
            y = fft(y,[],2);
            for i = 1:Nt
                z(:,i) = -1i*h(i)*(C*y(:,i)); % Frequency domain residual
            end
            
            %% Frequency-Domain Part
            for i = 1:Nt
              	%y(:,i) = y(:,i) + M2{i} \ z(:,i); % Apply second preconditioner, y = x1+x2
                y(:,i) = y(:,i) + (M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ z(:,i))))));
                z(:,i) = 1i*h(i)*(C*y(:,i)); % Frequency domain part of MVP
            end
            y = ifft(y,[],2,'symmetric');
            z = ifft(z,[],2,'symmetric');
            
            for i = 1:Nt
                j = 2*i - 1;
                x(:,i) = z(:,i) + (K{j}*y(:,i)) + (G{i}*y(:,i));
            end
            
            x = reshape(x,Nx*Nt,1);
        end
        
        function x = PC0_BDTF(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            z = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                %y(:,i) = M1{i} \ x(:,i); % Apply first preconditioner, y <- x1
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
            end
            
            y = fft(y,[],2) / Nt;
            for i = 1:Nt
                z(:,i) = -1i*h(i)*(C*y(:,i)); % Frequency domain residual
            end
            
            %% Frequency-Domain Part
            for i = 1:Nt
                %x(:,i) = y(:,i) + M2{i} \ z(:,i); % Apply second preconditioner, y = x1+x2
                x(:,i) = y(:,i) + (M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ z(:,i))))));
            end
        end
        
      	function x = MVP_PC0_BDTF_SUPER(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            z = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                j = 2*i - 1;
                %y(:,i) = M1{i} \ x(:,i); % Apply first preconditioner, y <- x1
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
                y(:,i) = (K{j}' * y(:,i)) + (G{i}' * y(:,i));
                z(:,i) = x(:,i) - (K{j} * y(:,i)) - (G{i} * y(:,i)); % Time domain part of residual
            end
            
            y = fft(y,[],2);
            z = fft(z,[],2);
            for i = 1:Nt
                z(:,i) = z(:,i)-1i*h(i)*(C*y(:,i)); % Frequency domain part of residual
            end
            
            %% Frequency-Domain Part
            for i = 1:Nt
                %y(:,i) = y(:,i) + (M2{i} \ z(:,i)); % Apply second preconditioner, z = x2
                y(:,i) = y(:,i) + (M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ z(:,i))))));
                z(:,i) = 1i*h(i)*(C*y(:,i)); % Frequency domain part of MVP
            end
            y = ifft(y,[],2,'symmetric');
            z = ifft(z,[],2,'symmetric');
            
            for i = 1:Nt
                j = 2*i - 1;
                x(:,i) = z(:,i) + (K{j}*y(:,i)) + (G{i}*y(:,i));
            end
            
            x = reshape(x,Nx*Nt,1);
        end
        
        function x = PC0_BDTF_SUPER(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            z = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                j = 2*i - 1;
                %y(:,i) = M1{i} \ x(:,i); % Apply first preconditioner, y <- x1
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
                y(:,i) = (K{j}' * y(:,i)) + (G{i}' * y(:,i));
                z(:,i) = x(:,i) - (K{j} * y(:,i)) - (G{i} * y(:,i)); % Time domain part of residual
            end
            
            y = fft(y,[],2);
            z = fft(z,[],2);
            for i = 1:Nt
                z(:,i) = z(:,i)-1i*h(i)*(C*y(:,i)); % Frequency domain part of residual
            end
            
            %% Frequency-Domain Part
            for i = 1:Nt
                %x(:,i) = y(:,i) + (M2{i} \ z(:,i));
                x(:,i) = y(:,i) + (M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ z(:,i))))));
            end
            x = x / Nt;
        end
        
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = HarmonicBalance(varargin{:});
            end
        end
    end
end