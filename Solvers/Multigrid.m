classdef IterativeHarmonicBalance < Solver
    properties
        Tolerance      = Static.setProperty(sqrt(eps));
        MaxItterations = Static.setProperty(1000);
    end
    
    methods
        function this = IterativeHarmonicBalance(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        function this = set.Tolerance(this,tol)
            this.Tolerance = this.setProperty(tol);
        end  
        
        function this = set.MaxItterations(this,maxItt)
            this.MaxItterations = this.setProperty(maxItt);
        end
        
        function solution = solve(this, model, x0)
            %% configure matrices
            matrixFactory = DynamicMatrixFactory(model);
            this.Matrices = matrixFactory;
            
            %% configure algorithm
          	maxItt       = this.MaxItterations.Value;
            errTol       = this.Tolerance.Value;
            shouldReport = this.ReportProgress;
            
            %% get times points
            t = matrixFactory.getTimePoints(this.TimePoints.Value);
            if this.TimePoints.Value == 1
                t  = t(1);
                Nt = 1;
            else
                %Levels = ceil(log2(this.TimePoints.Value / 8)) + 1;
                Levels = 6;
                Coarse = 1;
                Nt     = Coarse * 2 ^ (Levels - 1) + 1;
                t      = linspace(0, t(end), Nt);
            end

            this.Times = t;
            x          = cell(1, Nt);

            %% for each time point
            tic
            newtonError = inf;
            for i = 1:Nt
                %% create constant matrices
                f      = matrixFactory.f(t(i), 1);
                K      = matrixFactory.K(t(i), 1);

                %% get initial guess
                if i > 1
                    x{i} = x{i-1};
                else
                    if nargin == 3
                        x{1} = x0;
                    else
                        x{i} = zeros(size(f));
                    end
                end
                
                %% perform newton - raphson itteration
                jItt   = 1;
                relErr = 1;
                
                if shouldReport
                    TimeStepNumber = i
                end
                
                while jItt < maxItt && relErr > errTol
                    [G, g] = matrixFactory.G(t(i), x{i});
                    
                    r      = (K*x{i} + g - f);
                    J      = K + G;
                    
                    dx     = J \ r;
                    
                    x{i}   = (x{i} - dx);
                    
                    jItt   = jItt + 1;
                    relErr = norm(r) / norm(g - f);
                end
            end
            
            StaticAnalysisTime = toc
            
            %% evaluate residual
            Nt = Nt - 1;
            x  = x(1:Nt);
            r  = cell(size(x));
            f  = cell(size(x));
            g  = cell(size(x));
            G  = cell(size(x));
            K  = cell(size(x));
            
            J  = cell(1, Levels);           	
            C  = cell(1, Levels);
            
            h         = t(end) / Nt / 100;
            C{Levels} = cell(1, Nt);
            J{Levels} = cell(1, Nt);
            for i = 1:Nt
                k = 1;
                f{i}         = matrixFactory.f(t(i), h);
                K{i}         = matrixFactory.K(t(k), h);
                C{Levels}{i} = matrixFactory.C(t(k), h);
            end

            newtonError = inf;
            while newtonError > errTol
                for i = 1:Nt
                    k = 1;
                    [G{i}, g{i}] = matrixFactory.G(t(k), x{k});
                    j            = mod(i - 2, Nt) + 1;
                    r{i}         = f{i} - (K{i} * x{i} + C{Levels}{i} * (x{i} - x{j}) + g{i});
                    J{Levels}{i} = K{k} + G{k} + C{Levels}{k};
                    
%                     if i == 1
%                         S = sum(abs(J{Levels}{1}));
%                         S = 1 ./ sqrt(S);
%                         S = sparse(1:numel(S), 1:numel(S), S);
%                         [~,~,~,S] = ldl(J{Levels}{1});
%                         S = speye(size(K{1}));
%                     end
                    
                    J{Levels}{i} = J{Levels}{i};
                    C{Levels}{i} = C{Levels}{i};
                end
                
                %% Galerkin-Coarsening
%                 for i = (Levels - 1):-1:1
%                     n    = Nt / 2^(Levels - i);
%                     J{i} = cell(1, n);
%                     C{i} = cell(1, n);
%                     for j = 1:n
%                         l       = mod(2 * j - 2, 2 * n) + 1;
%                         k       = mod(2 * j - 3, 2 * n) + 1;
%                         J{i}{j} = J{i+1}{l} / 2 + J{i+1}{k} / 2;
%                         C{i}{j} = C{i+1}{k} / 2;
% 
% %                         J{i}{j} = J{i+1}{l} / 2 - C{i+1}{l} / 4 + J{i+1}{k} / 4;
% %                         
% %                         C{i}{j} = C{i+1}{k} / 2 - J{i+1}{k} / 4 + C{i+1}{l} / 4;
%                     end
%                 end
                
                %% Standard-Coarsening              
                for i = (Levels - 1):-1:1
                    n    = Nt / 2^(Levels - i);
                    J{i} = cell(1, n);
                    C{i} = cell(1, n);
                    for j = 1:n
                        k       = Nt / n * j - 1;
                        k       = 1;
                        h_ij    = h * Nt / n;
                        G_ij    = matrixFactory.G(t(k), x{k});
                        K_ij    = matrixFactory.K(t(k), h_ij);
%                         for l = 1:(2^(Levels - i) - 1)
%                             m       = mod(k - l, Nt) + 1;
%                             G_ij    = G_ij + matrixFactory.G(t(m), x{m});
%                             K_ij    = K_ij + matrixFactory.K(t(m), h_ij);
%                         end
%                         G_ij    = G_ij * 2^(i - Levels);
%                         K_ij    = K_ij * 2^(i - Levels);
                        C_ij    = matrixFactory.C(t(k), h_ij);
                        J{i}{j} = (K_ij + C_ij + G_ij);
                        C{i}{j} = C_ij;
                    end
                end
                
                tic
                eigMax = zeros(1, Levels);
                optRP  = zeros(1, Levels);
               	L      = cell(1, Levels);
                D      = cell(1, Levels);
               	P      = cell(1, Levels);
                S      = cell(1, Levels);
                for i = 1:numel(J)
                    for j = 1:numel(J{i})
                        [L{i}{j}, D{i}{j}, P{i}{j}, S{i}{j}] = ldl(J{i}{j});
                    end
                end
                FactorizationTime = toc
                
                
                tic
                for i = numel(J):-1:2
                    eigFun    = @(var)(S{i}{1} * (P{i}{1} * (L{i}{1}.' \ (D{i}{1} \ (L{i}{1} \ (P{i}{1}.' * (S{i}{1} * ( - C{i}{1} * var))))))));
%                     eigFun    = @(var)(- (C{i}{1} * S{i}{1} * (P{i}{1} * (L{i}{1}.' \ (D{i}{1} \ (L{i}{1} \ (P{i}{1}.' * (S{i}{1} * var))))))));
                    eigFun    = @(var)(J{end}{1}\(C{end}{1}.'*(C{end}{1}*(J{end}{1} \ var))));
                    eigMax(i) = eigs(eigFun, length(L{i}{1}), 1, 'LM');
                    optRP(i)  = 1 / (1 + eigMax(i).^2);
                end
                RelaxationParameterTime = toc
                
                r = cell2mat(r);
                f = cell2mat(f);
                g = cell2mat(g);

                newtonError = norm(reshape(r, [], 1)) / norm(reshape(g, [], 1) - reshape(f, [], 1))
                itterError  = newtonError;
                
                %% Matrix-Vector Product
                Afun = @(var)(mvp(this, var, J{Levels}, C{Levels}));
                
                %% Multigrid Preconditioner
                Mfun = @(var)(mgcyc(this, Levels, [], 0 * var , L, D, P, S, C, J, var, 2, 2, optRP));
                
                %% Block-Diagonal Preconditioner
%                 Mfun = @(x)(bdiagpc(this, x, L{Levels}, U{Levels}, P{Levels}, Q{Levels}, R{Levels}));
                
                Afun = @(var)(Afun(Mfun(var)));
                
                tic
                r = 1 - 2 * rand(numel(r),1);
                r = r / norm(r);
                dx = gmres(Afun, r, 2, [], 1);%, Mfun);
                norm(Afun(dx)-r) / norm(r)
                GMRESTime = toc
                
                temp = [J{2}{1} -C{2}{2};-C{2}{1} J{2}{2}];
                
                
                dx = reshape(dx, [], Nt);
                
%                 dx = 0 * r;
%                 dy = 0 * r;
%                 s  = 0 * r;
%                 while itterError > errTol
%                     for i = 1:Nt
%                         j = mod(i - 2, Nt) + 1;
%                         dy(:, i) = P{i}.'\(L{i}.'\(D{i}\(L{i}\(P{i}\(C * dx(:,j) + r(:,i))))));
%                     end
%                     
%                     for i = 1:Nt
%                         j       = mod(i - 2, Nt) + 1;
%                         s(:, i) = J{i} * dy(:, i) - C * dy(:, j) - r(:, i);
%                     end
%                     
%                     itterError = norm(reshape(s,[],1)) / norm(reshape(r,[],1))
%                     
%                     dx = dy;
%                 end
                
              	r = mat2cell(r, length(r), ones(1, Nt));
                f = mat2cell(f, length(f), ones(1, Nt));
                g = mat2cell(g, length(g), ones(1, Nt));
            end
            
            x(end+1) = x(1);
            Nt       = Nt + 1;
            x_t      = x;
            for i = 1:Nt
                j = mod(i - 2, Nt) + 1;
                x_t{i} = (x{i} - x{j}) / h;
            end

            if shouldReport
                toc;
            end

            %% Save Solution
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x, x_t);
            solution           = Solution(this);
        end
        
        function u = mgcyc(this, k, gamma, u, L, D, P, S, C, J, f, upsilon1, upsilon2, optRP)
            if norm(f) == 0
                u = f;
            else
                N = 2;
                if k > N+1
                    u = smooth(this, upsilon1, u, L{k}, D{k}, P{k}, S{k}, C{k}, f, optRP(k), J{k});
                    d = defect(this, f, J{k}, C{k}, u);
                    d = restrict(this, d, numel(J{k}));
                    v = mgcyc(this, k-1, gamma, 0 * d, L, D, P, S, C, J, d, upsilon1, upsilon2, optRP);
                    v = interpolate(this, v, numel(J{k}));
                    u = u + v;
                    u = smooth(this, upsilon2, u, L{k}, D{k}, P{k}, S{k}, C{k}, f, optRP(k), J{k});
                else
                    temp      = cell(2^N, 2^N);
                    [temp{:}] = deal(sparse(length(J{k}{1}), length(J{k}{1})));
                    for i = 1:(2^N)
                        temp{i,i} = J{k}{i};
                        j         = mod(i-2, 2^N) + 1;
                        temp{i,j} = - C{k}{i};
                    end
                    u = cell2mat(temp) \ reshape(f,[],1);
                end
            end
            
            N = 8;
            M = abs(K{1}-K{1}.') > max(max(abs(K{1}))) * sqrt(eps);
            K = (sum(M,2) > 1) | (diag(K{1}) == 0);
            I = ~K;
            
            
            JinvC = J{3}{1}(I,I) \ C{3}{1}(I,I);
            II    = speye(size(JinvC));
            
            w      = 1;
            Z      = cell(4,4);
            [Z{:}] = deal(sparse(length(II),length(II)));
            Z{1,1} = (1-w)*II;
            Z{1,3} = (w*JinvC)^2;
            Z{1,4} = (1-w)*w*JinvC;
            Z{2,1} = w*JinvC;
            Z{2,2} = (1-w)*II;
            Z{3,1} = (w*JinvC)^2;
            Z{3,2} = (1-w)*w*JinvC;
            Z{3,3} = (1-w)*II;
            Z{4,3} = w*JinvC;
            Z{4,4} = (1-w)*II;
            
            Z = cell2mat(Z);
            lz = eigs(Z,1000,'LM');
            
            w      = 0.5;
           	Q      = cell(4,4);
            [Q{:}] = deal(sparse(length(II),length(II)));
            Q{1,1} = (1-w)*II;
            Q{1,4} = w*JinvC;
            Q{2,2} = (1-w)*II;
            Q{2,1} = w*JinvC;
            Q{3,3} = (1-w)*II;
            Q{3,2} = w*JinvC;
            Q{4,4} = (1-w)*II;
            Q{4,3} = w*JinvC;
            
            Q = cell2mat(Q);
            
%             Z = cell(2*N,2*N);
%             Ahat = cell(N,N);
%             [Ahat{:}] = deal(sparse(sum(I),sum(I)));
%             n = [sum(I)*ones(1,N),sum(K)*ones(1,N)];
%             for i = 1:(2*N)
%                 for j = 1:(2*N)
%                     Z{i,j} = sparse(n(i), n(j));
%                 end
%             end
%             
%             for i = 1:N
%                 j          = mod(i-2,N) + 1;
%                 Ahat{i,i}  = J{log2(N)+1}{1}(I,I);
%                 Z{i, i}    = J{log2(N)+1}{1}(I,I);
%                 Z{i, j}    = -C{log2(N)+1}{1}(I,I);
%                 Z{i, i+N}  = J{log2(N)+1}{1}(I,K);
%                 Z{i, j+N}  = -C{log2(N)+1}{1}(I,K);
%                 Z{i+N,i}   = J{log2(N)+1}{1}(K,I);
%                 Z{i+N,j}   = -C{log2(N)+1}{1}(K,I);
%                 Z{i+N,i+N} = J{log2(N)+1}{1}(K,K);
%                 Z{i+N,j+N} = -C{log2(N)+1}{1}(K,K);
%             end
%             Ahat = cell2mat(Ahat);
%             A = cell2mat(Z(1:N,1:N));
%             U = cell2mat(Z(1:N,(N+1):(2*N)));
%             L = cell2mat(Z((N+1):(2*N),1:N));
%             B = cell2mat(Z((N+1):(2*N),(N+1):(2*N)));
%             A2 = A-U*(B\L);
%             
%             T = [(U*(B\L)-A)/Ahat2, -((speye(size(A))-A/Ahat2)*U/B + U*(B\L)*(Ahat2\U)/B);
%                  sparse(length(B),length(A)), sparse(length(B),length(B))];
        end
        
        function u = smooth(this, upsilon, u, L, D, P, S, C, f, optRP, J)
            Nt = numel(L);
            u  = reshape(u, [], Nt);
            f  = reshape(f, [], Nt);
            
            for k = 1:upsilon
                d = defect(this, f, J, C, u);
                d = reshape(d, [], Nt);
                for i = 1:Nt
                    u(:, i) = u(:, i) + optRP * (S{i} * (P{i} * (L{i}.' \ (D{i} \ (L{i} \ (P{i}.' * (S{i} * (d(:, i)))))))));
                end
            end
            
            u = reshape(u, [], 1);
        end
        
        function d = defect(this, f, J, C, u)
            Nt = numel(J);
            f  = reshape(f, [], Nt);
            u  = reshape(u, [], Nt);
            d  = u;
            for i = 1:Nt
                j       = mod(i - 2, Nt) + 1;
                d(:, i) = f(:, i) - (J{i} * u(:, i) - C{j} * u(:, j));
            end
            d  = reshape(d, [], 1);
        end
        
        function d_c = restrict(this, d_f, N)
            %% Injection
            d_f = reshape(d_f, [], N);
            d_c = d_f(:,1:2:end);
            d_c = reshape(d_c, [], 1);

            %% Backward-Euler Restriction
%             d_f = reshape(d_f, [], N);
%             
%             I   = 1:2:(N-1);
%             d_c = d_f(:, I);
%             
%             I   = [N, 2:2:(N-2)];
%             d_c = d_c + d_f(:, I);
%             
%             d_c = d_c / 2;
%             d_c = reshape(d_c, [], 1);

            %% Full Weighting
%             d_f = reshape(d_f, [], N);
%             d_c = 2 * d_f(:,1:2:end);
%             d_c = d_c + d_f(:,[end, 2:2:(end-2)]);
%             d_c = d_c + d_f(:,[2:2:(end-2), 1]);
%             d_c = d_c / 4;
%             d_c = reshape(d_c, [], 1);

            %% DFT
%             d_f            = reshape(d_f, [], N);
%             if N > 2
%                 d_f            = fft(d_f, [], 2) / N;
%                 d_c            = d_f(:,[1:(N/4) (3*N/4+1):N]);
%                 d_c(:,N/4 + 1) = d_c(:,N/4 + 1) * 0;
%                 d_c            = ifft(d_c, [], 2) * N / 2;
%                 d_c            = reshape(d_c, [], 1);
%             else
%                 d_c = d_f(:, 1) / 2 + d_f(:, 2) / 2;
%             end
        end
        
        function v_f = interpolate(this, v_c, N)
            %% Linear Interpolation
            v_c       = reshape(v_c, [], N / 2);
            v_f       = zeros(size(v_c) * diag([1 2]));
            
            I         = 1:2:(N - 1);
            v_f(:, I) = v_c;
            
            I         = I + 1;
            v_f(:, I) = (v_c + v_c(:,[2:end 1])) / 2; 
            
            v_f       = reshape(v_f, [], 1);

            %% DFT
%             if N > 2
%                 v_c                 = reshape(v_c, [], N / 2);
%                 v_f                 = zeros(size(v_c) * diag([1 2])); 
%                 v_c                 = fft(v_c, [], 2) * 2 / N;
%                 v_f(:, 1:(N/4))     = v_c(:, 1:(N/4));
%                 v_f(:, (3*N/4+2):N) = v_c(:, (N/4+2):(N/2));
%                 v_f(:, N / 2 + 1)   = v_c(:, N / 4 + 1);
%                 v_f                 = ifft(v_f, [], 2) * N;
%                 v_f                 = reshape(v_f, [], 1);
%             elseif N == 2
%                 v_f = [v_c; v_c];
%             end 
        end
        
        function y = mvp(this, x, J, C)
            Nt = numel(J);
            x  = reshape(x, [], Nt);
            y  = x;
            for i = 1:Nt
                j       = mod(i - 2, Nt) + 1;
                y(:, i) = J{i} * x(:, i) - C{j} * x(:, j);
            end
            y  = reshape(y, [], 1);
        end
        
        function x = bdiagpc(this, b, L, D, P, S)
            Nt = numel(L);
            b  = reshape(b, [], Nt);
            x = 0 * b;
            for i = 1:Nt
                x(:, i) = S{i} * (P{i} * (L{i}.' \ (D{i} \ (L{i} \ (P{i}.' * (S{i} * (b(:, i))))))));
            end
            x = reshape(x, [], 1);
        end
    end
    
    methods (Static)        
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = IterativeHarmonicBalance(varargin{:});
            end
        end
    end
end