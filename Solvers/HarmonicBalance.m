classdef HarmonicBalance < Solver
    properties
        AdaptiveTol   = 1e-2;
        GMRESTol      = 1e-2;
        NewtonTol     = 1e-1;
        ColocationTol = 1e-4;
        
        MinAdaptiveIter = 0;
        MaxAdaptiveIter = [];
        MaxGMRESIter   	= 100;
        MinNewtonIter  	= 1;
        MaxNewtonIter  	= 20;
        
        StorageLevel = 3;
        
        Strategy = 'target';
        Plan = [];
    end
    
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
        function solution = solve(this, model, x0)
            %% Setup Matrices
            getMatrix = DynamicMatrixFactory(model);
            this.Matrices = getMatrix;
            Nx = length(getMatrix.f(0));

            t = model.getTimePoints(this.TimePoints);
            
            %% Initialize
            switch lower(this.Strategy)
                case 'target'
                    if nargin < 3
                        plan = factor((numel(t) - 1)/6);
                        plan = sort(plan, 'descend');
                        plan = [3, 2, plan, 1];
                    else
                        plan = factor((numel(t) - 1) / size(x0,2));
                        plan = sort(plan, 'descend');
                        plan = [size(y0,2), plan, 1];
                    end
                case 'plan'
                    plan = this.Plan;
                    if nargin == 3 && plan(1) ~= size(x0,2)
                        plan = [size(x0,2), plan];
                    end
                    if plan(end) ~= 1
                        plan(end+1) = 1;
                    end
                otherwise
                    error('MotorProtor:HarmonicBalance','Unknown strategy %s. Try "target" or "plan."',this.Strategy);
            end
            this.Plan = plan;
            
            Nt = plan(1);
            Nh = ceil(Nt/2)-1;
            
            T = t(end);
            t = linspace(0, T, Nt+1);
            
            adaptTol  = this.AdaptiveTol;
            gmresTol  = this.GMRESTol;
            
            gmresIter = this.MaxGMRESIter;
            
            minAIter = this.MinAdaptiveIter;
            if isempty(this.MaxAdaptiveIter)
                this.MaxAdaptiveIter = length(plan);
            end
            maxAIter = this.MaxAdaptiveIter;
            
            store   = this.StorageLevel;
            verbose = this.Verbose;
            
            %% Allocate
            tic
            
            if nargin < 3
                X = zeros(Nx, Nt);
            else
                X = fft(x0,[],2) / Nt;
            end
                
            C = getMatrix.C(0,1,1);
            
            K = cell(1,Nt);
            f = zeros(Nx,Nt);
            
            for i = 1:Nt
                K{i} = getMatrix.K(t(i),1);
                f(:,i) = getMatrix.f(t(i),1);
            end
            
            switch store
                case 3
                    MVP = @HarmonicBalance.MVPStoredLU;
                    PC = @HarmonicBalance.PCStoredLU;
                case 2
                    MVP = @HarmonicBalance.MVP2;
                    PC = @HarmonicBalance.PC2;
                otherwise
                    error('MotorProto:HarmonicBalance','No Implementation');
            end
            
            %% Start Solve
            aIter  = 1;
            vareps = inf;
            final  = false;
            while ~final || ((vareps > adaptTol) && ~final)
                %% Subdivide Time Interval
                if ((aIter < maxAIter) && (vareps > adaptTol)) || (aIter < minAIter)
                    aIter = aIter + 1;
                    final = final || (aIter == length(plan));
                   	d = plan(aIter);
                    
                    if mod(Nt,2) == 0 && d > 1
                        X = [X(:,1:Nh), X(:,Nh+1) / 2, zeros(Nx, (d-1)*Nt-1), conj(X(:,Nh+1)) / 2, X(:,(Nh+2):Nt)];
                    else
                        X = [X(:,1:(Nh+1)), zeros(Nx, (d-1)*Nt), X(:,(Nh+2):Nt)];
                    end
                    
                    Nt = d*Nt;
                    Nh = ceil(Nt/2)-1;
                    if mod(Nt,2) == 0
                        h = [0:Nh, 0, -Nh:-1] * 2 * pi / T;
                    else
                        h = [0:Nh, -Nh:-1] * 2 * pi / T;
                    end
                    
                   	w = (h.^(6) ./ (h.^6+(pi*Nt/2/T)^6)) .^(1/2);
                    if mod(Nt,2) == 0
                        w(Nt/2) = 1;
                    end
                    
                    if verbose
                        display(sprintf('Time Points = %d\n',Nt));
                    end
                    
                    t = linspace(0,T,Nt+1);
                    
                    G = cell(1,Nt);
                    g = zeros(Nx,Nt);
                    
                    K = [K, cell(1, Nt*(d-1))];
                    f = [f, zeros(Nx, Nt*(d-1))];
                    
                    for i = (Nt-d+1):-d:1
                        j      = (i+d-1) / d;
                        K{i}   = K{j};
                        f(:,i) = f(:,j);
                    end
                    
                    for i = setdiff(1:Nt, 1:d:Nt)
                        K{i}   = getMatrix.K(t(i),1);
                        f(:,i) = getMatrix.f(t(i),1);
                    end
                    
                    R = zeros(Nx,Nt);
                    
                    if store == 3
                        M1 = cell(5,Nt);
                        M2 = cell(5,Nt);
                    else
                        M1 = cell(1,Nt);
                        M2 = cell(1,Nt);
                    end
                    
                    nIter = 0;
                else
                    final = true;
                end
                
                if verbose && final
                    display(sprintf('Colocating to a tolerance of %f with %d time points.\n', min(this.AdaptiveTol, this.ColocationTol), Nt));
                end
                    
                %% Calculate Residual Vector
                vareps = 0;
                alpha  = 0;
                beta   = 0;
                
                for i = 1:Nt
                    R(:,i) = 1i*h(i)*(C*X(:,i));
                    alpha  = alpha + norm(X(:,i))^2;
                    vareps = vareps + w(i)^2*norm(X(:,i))^2;
                end
                if alpha ~= 0
                    vareps = vareps / alpha;
                    alpha  = 0;
                end
                
                r = ifft(R,[],2,'symmetric') * Nt;
                x = ifft(X,[],2,'symmetric') * Nt;
                for i = 1:Nt
                    [G{i},g(:,i)] = getMatrix.G(t(i), x(:,i), 1);
                    r(:,i) = r(:,i) + K{i} * x(:,i) + g(:,i);
                    
                    beta = beta + 0.5*norm(r(:,i)+f(:,i))^2;
                    
                    r(:,i) = r(:,i) - f(:,i);
                    
                    alpha = alpha + norm(r(:,i))^2;
                end
                vareps = sqrt(vareps);
                alpha  = sqrt(alpha);
                beta   = sqrt(beta);
                
                if (vareps < this.AdaptiveTol) || final
                	tol = min(this.AdaptiveTol, this.ColocationTol);
                else
                    eta = max(0,(log10(this.AdaptiveTol)-log10(vareps)) / log10(this.AdaptiveTol));
                    tol = 10^(log10(vareps) * eta + log10(this.ColocationTol) * (1-eta));
                  	tol = max(tol, vareps*this.NewtonTol);
                    nIter = 0;
                end
                
                %% Start Newton Iteration
                while ((alpha > tol * beta) && (nIter < this.MaxNewtonIter)) || (nIter < this.MinNewtonIter)
                    %% Solve Linearized Equation
                    %% Average of Static Part
                    J0 = K{1} + G{1};
                    for i = 2:Nt
                        J0 = J0 + K{i} + G{i};
                    end
                    J0 = J0 / Nt;
                    
                    %% Calculate and Factor both Preconditioners
                    for i = 1:Nt
                        M1{1,i} = K{i}+G{i}+(Nt/T)*C;
                        M2{1,i} = J0 + 1i*h(i)*C;
                        
                        if store == 3
                            [M1{3,i},M1{4,i},M1{2,i},M1{5,i},M1{1,i}] = lu(M1{1,i});
                            [M2{3,i},M2{4,i},M2{2,i},M2{5,i},M2{1,i}] = lu(M2{1,i});
                        end
                    end
                    
                    %% GMRES Solve
                    dX = gmres(MVP, r(:), gmresIter, gmresTol, 1, [], [], r(:), K, G, C, M1, M2, h);
                    dX = PC(dX,K,G,C,M1,M2,h);
                    X  = X - dX;
                    
                    %% Calculate Residual Vector
                    vareps = 0;
                    alpha  = 0;
                    beta   = 0;
                    
                    for i = 1:Nt
                        R(:,i) = 1i*h(i)*(C*X(:,i));
                        alpha  = alpha + norm(X(:,i))^2;
                        vareps = vareps + w(i)^2*norm(X(:,i))^2;
                    end
                    
                    if alpha ~= 0
                        vareps = vareps / alpha;
                        alpha  = 0;
                    end
                    
                    r = ifft(R,[],2,'symmetric') * Nt;
                    x = ifft(X,[],2,'symmetric') * Nt;
                    for i = 1:Nt
                        [G{i},g(:,i)] = getMatrix.G(t(i), x(:,i), 1);
                        r(:,i) = r(:,i) + K{i} * x(:,i) + g(:,i);
                        
                        beta = beta + 0.5*norm(r(:,i)+f(:,i))^2;
                        
                        r(:,i) = r(:,i) - f(:,i);
                        
                        alpha = alpha + norm(r(:,i))^2;
                    end
                    vareps = sqrt(vareps);
                    alpha  = sqrt(alpha);
                    beta   = sqrt(beta);
                    
                    if (vareps < this.AdaptiveTol) || final
                        tol = min(this.AdaptiveTol, this.ColocationTol);
                    else
                        eta = max(0,(log10(this.AdaptiveTol)-log10(vareps)) / log10(this.AdaptiveTol));
                        tol = 10^(log10(vareps) * eta + log10(this.ColocationTol) * (1-eta));
                        tol = max(tol, vareps*this.NewtonTol);
                    end
                    
                    nIter = nIter + 1;
                    
                    if verbose
                        display(sprintf('Iteration %d, Discrete Residual = %f, Tolerance = %f, Discretization Error = %f, Tolerance = %f \n', nIter, alpha / beta, tol, vareps, adaptTol));
                    end
                end
            end
            this.SimulationTime = toc;
            
           	if this.Verbose
                display(sprintf('Simulation Time = %0.3g seconds\n', this.SimulationTime));
            end
            
            %% Post Processing
            x = ifft(X, [], 2, 'symmetric') * Nt;
            x_t = X;
            for i = 1:Nt
                x_t(:,i) = 1i*h(i)*x_t(:,i);
            end
            x_t = ifft(x_t, [], 2, 'symmetric') * Nt;
            
            x = [x, x(:,1)];
            x_t = [x_t, x_t(:,1)];
            
            x = mat2cell(x, Nx, ones(1,Nt+1));
            x_t = mat2cell(x_t, Nx, ones(1,Nt+1));
            
            [x, x_t] = getMatrix.doPostProcessing(x, x_t);
            this.Times = linspace(0,T,Nt+1);
            this.Y   = X;
            this.X   = x;
            this.X_t = x_t;
            
            solution = Solution(this);
        end
    end
    
    methods (Static)
    	function x = MVPStoredLU(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
                x(:,i) = x(:,i) - (K{i}*y(:,i)) - (G{i}*y(:,i));
            end
            
            %% Frequency-Domain Part
            x = fft(x,[],2);
            y = fft(y,[],2);
            for i = 1:Nt
                x(:,i) = x(:,i)-1i*h(i)*(C*y(:,i)); % Frequency domain residual
                y(:,i) = y(:,i)+(M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ x(:,i)))))); % Apply second preconditioner, y = x1+x2
                x(:,i) = 1i*h(i)*(C*y(:,i)); % Frequency domain part of MVP
            end
            
            %% Time-Domain Part
            x = ifft(x,[],2,'symmetric');
            y = ifft(y,[],2,'symmetric');
            for i = 1:Nt
                x(:,i) = x(:,i) + (K{i}*y(:,i)) + (G{i}*y(:,i));
            end
            
            x = reshape(x,Nx*Nt,1);
        end
        
        function y = PCStoredLU(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                y(:,i) = M1{5,i} * (M1{4,i} \ (M1{3,i} \ (M1{2,i} * (M1{1,i} \ x(:,i)))));
                x(:,i) = x(:,i) - (K{i}*y(:,i)) - (G{i}*y(:,i));
            end
            
            %% Frequency-Domain Part
            x = fft(x,[],2) / Nt;
            y = fft(y,[],2) / Nt;
            for i = 1:Nt
                x(:,i) = x(:,i)-1i*h(i)*(C*y(:,i)); % Frequency domain residual
                y(:,i) = y(:,i)+(M2{5,i} * (M2{4,i} \ (M2{3,i} \ (M2{2,i} * (M2{1,i} \ x(:,i)))))); % Apply second preconditioner, y = x1+x2
            end
        end
        
     	function x = MVP2(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);
            
            %% Time-Domain Part
            for i = 1:Nt
                y(:,i) = M1{i}\x(:,i);
                x(:,i) = x(:,i) - (K{i}*y(:,i)) - (G{i}*y(:,i));
            end
            
            %% Frequency-Domain Part
            x = fft(x,[],2);
            y = fft(y,[],2);
            for i = 1:Nt
                x(:,i) = x(:,i) - 1i*h(i)*(C*y(:,i)); % Frequency domain residual
                y(:,i) = y(:,i) + M2{i}\x(:,i); % Apply second preconditioner, y = x1+x2
                x(:,i) = 1i*h(i)*(C*y(:,i)); % Frequency domain part of MVP
            end
            
            %% Time-Domain Part
            x = ifft(x,[],2,'symmetric');
            y = ifft(y,[],2,'symmetric');
            for i = 1:Nt
                x(:,i) = x(:,i) + (K{i}*y(:,i)) + (G{i}*y(:,i));
            end
            
            x = reshape(x,Nx*Nt,1);
        end
        
        function y = PC2(x,K,G,C,M1,M2,h)
            Nt = length(G);
            Nx = length(x) / Nt;
            
            x = reshape(x,Nx,Nt);
            y = zeros(Nx,Nt);f
            
            %% Time-Domain Parts
            for i = 1:Nt
                y(:,i) = M1{i}\x(:,i);
                x(:,i) = x(:,i) - (K{i}*y(:,i)) - (G{i}*y(:,i));
            end
            
            %% Frequency-Domain Part
            x = fft(x,[],2) / Nt;
            y = fft(y,[],2) / Nt;
            for i = 1:Nt
                x(:,i) = x(:,i) - 1i*h(i)*(C*y(:,i)); % Frequency domain residual
                y(:,i) = y(:,i) + M2{i}\x(:,i); % Apply second preconditioner, y = x1+x2
            end
        end
        
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = HarmonicBalance(varargin{:});
            end
        end
    end
end