classdef HarmonicBalanceDomainDecomposition < Solver
    %HarmonicBalanceDD.m Harmonic balance periodic steady-state solver based on a domain decomposition algorithm
    %
    % HarmonicBalanceDD properties:
    %   NewtonTolerance   - Tolerance of outer Newton-Raphson iterations
    %   NewtonIterations  - Maximum number of Newton-Raphson iterations
    %   GMRESTolerance    - Tolerance of inner GMRES iterations
    %   Adaptive          - Use adaptive refinement
    %   DirectSolver      - Use directive solver
    %
    % HarmonicBalanceDD inherits properties and methods Solver.
    %
    % See also MotorProto, Solver
    
%{
properties:
 	%NewtonTolerance - Sets the relative residual tolerance
 	%	NewtonTolerance sets the relative residual tolerance for the Newton-Raphson
    %   iteration that occurs at each time step. The default value is sqrt(eps).
    %
    % See also Static, NewtonIterations
    NewtonTolerance;
%}
    
    properties
        NewtonTolerance   = 1e-6;
        NewtonIterations  = 100;
        GMRESTolerance    = 1e-6;
        Adaptive          = false;
        TimeStepTolerance = 0;
        DirectSolver      = false;
    end
    
    methods
        %% Constructor
        function this = HarmonicBalanceDomainDecomposition(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Solve
        function solution = solve(this, model,~)
            directSolver = this.DirectSolver;
            if ~directSolver
                %% Domain Decomposition
                if this.Adaptive
                    Nt        = 2.^(2:ceil(log2(this.TimePoints)));
                    nAdaptive = length(Nt);
                else
                    Nt        = this.TimePoints;
                    nAdaptive = 1;
                end
                
                %% Preprocessing
                matrixFactory = HarmonicBalanceDomainDecompositionMatrixFactory(copy(model));
                this.Matrices = matrixFactory;
                t  = model.getTimePoints(this.TimePoints);
                T  = t(end);
                f  = 1 / T;
                w  = 2 * pi * f;
                Nt = length(t);
                if mod(Nt,2) == 1
                    Nt = Nt - 1;
                end
                
                for iAdaptive = 1:nAdaptive
                    t      = linspace(0,T,Nt(iAdaptive)+1);
                    t(end) = [];
                    Nh     = Nt(iAdaptive) / 2;
                    omega  = 1i * w * [0:(Nh-1), 0, ((1-Nh):1:-1)];
                    Omega  = sparse(1:Nt(iAdaptive),1:Nt(iAdaptive),omega);

                    [K1,K2]    = matrixFactory.K(t);
                    C1         = matrixFactory.C;
                    f          = matrixFactory.f(t);

                       
                    K = this.addMatrices(K1,K2);

                    [I1,I2,Ib] = matrixFactory.getDomainVectors;

                    Ik      = I1 | I2;
                    
                    [Ek,~] = this.extensionMatrices(Ik,false & Ik, Ib);

                    [K1n,K1ii,K1ib,K1bi,K1bb] = this.splitMatrix(K1,I1,Ib);
                    [C1n,C1ii,C1ib,C1bi,C1bb] = this.splitMatrix(C1,I1,Ib);

                    [M1s] = this.factorFrequencyDomainMatrices(K1ii,C1ii,K1n,C1n,omega);

                    [Ko,Co,Po,Eo,E2To,To,Ro,T2Eo] = matrixFactory.OverlappingPreconditioner(t);

                    if iAdaptive == 1
                        x     = 0 * f;
                    else
                        y       = fft(x,[],2) / Nt(iAdaptive-1);
                        x       = zeros(length(y),Nt(iAdaptive));
                        Ix      = [1:1:(Nt(iAdaptive-1)/2), Nt(iAdaptive  )/2+1, fliplr(Nt(iAdaptive-1):-1:(Nt(iAdaptive-1)/2+2)) + Nt(iAdaptive-1)];
                        Iy      = [1:1:(Nt(iAdaptive-1)/2), Nt(iAdaptive-1)/2+1, fliplr(Nt(iAdaptive-1):-1:(Nt(iAdaptive-1)/2+2))];
                        x(:,Ix) = y(:,Iy);
                        x       = ifft(x,[],2) * Nt(iAdaptive);
                    end

                    [G,g] = matrixFactory.G(t,x);

                    r = this.calculateResidual(x,K1,K2,C1,Omega,g,f);

                    relres = norm(r(:)) / norm(f(:));
                    itter  = 0;
                    tic
                    while (itter < 20) && (relres > this.NewtonTolerance)
                        %% Factor Nonlinear
                        
                        Gk = this.addMatrices(G,K);
                        Go = this.averageMatrices(G);
                        G  = this.addMatrices(G,K2);

                        [K2n,K2ii,K2ib,K2bi,K2bb] = this.splitMatrix(G,I2,Ib);

                        [M2s] = this.factorTimeDomainMatrices(K2ii,K2n);

                        
                        [Gk,~,~,~,~] = this.splitMatrix(Gk,Ik,Ib);
                        
                        [~,Mk]       = this.factorTimeDomainMatrices([],Gk);
                        
                        Mo           = this.factorOverlapMatrices(Ko,Go,Co,Po,omega);

                        %% Lower Triangular Solve
                        [y1,y2,yb] = this.lowerTriangularSolve(r,I1,I2,Ib,K1ii,K1bi,C1ii,C1bi,omega,K2ii,K2bi);

                        %% GMRES Schur Complement Solve
                        A  = @(x)(this.MVP(x,M1s,K1ib,K1bi,K1bb,C1ii,C1ib,C1bi,C1bb,omega,M2s,K2ib,K2bi,K2bb));
                        M  = @(x)(this.PCk(x,Ek,Mk));
                        %M  = @(x)(0*x);
                        M  = @(x)(this.PCo(x,A,M,omega,Mo,Eo,E2To,To,Ro,T2Eo));
                        AM = @(x)(A(M(x)));
                        dy = gmres(AM,yb(:),[],this.GMRESTolerance,300);
                        dy = M(dy(:));
                        dy = real(dy);

                        %% Upper Triangular Solve
                        dx = this.upperTriangularSolve(y1,y2,dy,I1,I2,Ib,K1ii,K1ib,C1ii,C1ib,omega,K2ii,K2ib);

                        x     = x - dx;
                        [G,g] = matrixFactory.G(t,x);

                        r = this.calculateResidual(x,K1,K2,C1,Omega,g,f);

                        relres = norm(r(:)) / norm(f(:))
                        itter  = itter + 1;
                    end
                    this.Times = t;
                    this.SimulationTime = toc
                end
            elseif directSolver == 1
                t          = model.getTimePoints(this.TimePoints);
                t          = linspace(0,t(end),this.TimePoints+1);
                this.Times = t;

                T      = t(end);
                f      = 1 / T;
                w      = 2 * pi * f;
                t(end) = [];
                Nt     = numel(t);
                Nh     = Nt / 2;
                omega  = 1i * w * [0:(Nh-1), 0, ((1-Nh):1:-1)];
                Omega  = sparse(1:Nt,1:Nt,omega);
                
                %% Direct HB with DynamicMatrixFactory
                matrixFactory = DynamicMatrixFactory(copy(model));
                this.Matrices = matrixFactory;
                
                c = matrixFactory.C(0,1);
                K = cell(Nt,1);
                f = cell(1,Nt);
                for i = 1:Nt
                    K{i} = matrixFactory.K(t(i),1);
                    f{i} = matrixFactory.f(t(i),1);
                end
                K = blkdiag(K{:});
                f = [f{:}];
                
                CD  = cell(Nt,Nt);
                PCD = cell(Nt,Nt);
                [PCD{:}] = deal(sparse(size(c,1),size(c,2)));
                
                D = exp(-1i*2*pi*((0:(Nt-1)).'*(0:(Nt-1)))/Nt);
                P = D' / Nt;
                
                for i = 1:Nt
                    for j = 1:Nt
                        CD{i,j} = omega(i)*c*D(i,j);
                    end
                end
                
                for i = 1:Nt
                    for j = 1:Nt
                        for l = 1:Nt
                            PCD{i,j} = PCD{i,j} + P(i,l) * CD{l,j};
                        end
                    end
                end
                PCD = real(cell2mat(PCD));
                
                g = f * 0;
                f = f(:);
                x = f * 0;
                
                K  = K + PCD;
                Jm = cell(Nt,1);
                
                x = reshape(x,[],Nt);
                for i = 1:Nt
                   [Jm{i},g(:,i)] = matrixFactory.G(t(i),x(:,i));
                end
                x = x(:);
                r = K * x + g(:) - f;
                
                tic
                relres = norm(r) / norm(f)
                itter  = 0;
                while ((relres > sqrt(eps)) && (itter < 20)) || (itter == 0)
                    J  = K + blkdiag(Jm{:});
                    dy = J \ r;
                    x  = x - dy;
                    
                    x = reshape(x,[],Nt);
                    for i = 1:Nt
                       [Jm{i},g(:,i)] = matrixFactory.G(t(i),x(:,i));
                    end
                    x = x(:);
                    r = K * x + g(:) - f;
                    
                    relres = norm(r) / norm(f)
                    itter  = itter + 1
                end
                toc
                x = reshape(x,[],Nt);
            elseif directSolver == 2
                t          = model.getTimePoints(this.TimePoints);
                t          = linspace(0,t(end),this.TimePoints+1);
                this.Times = t;

                T      = t(end);
                f      = 1 / T;
                w      = 2 * pi * f;
                t(end) = [];
                Nt     = numel(t);
                Nh     = Nt / 2;
                omega  = 1i * w * [0:(Nh-1), 0, ((1-Nh):1:-1)];
                Omega  = sparse(1:Nt,1:Nt,omega);
                %% Direct HB with HarmonicBalanceDomainDecompositionMatrixFactory
                matrixFactory = HarmonicBalanceDomainDecompositionMatrixFactory(copy(model));
                this.Matrices = matrixFactory;
                
                c       = matrixFactory.C;
               	[k1,k2] = matrixFactory.K(t);
                K       = cell(1,Nt);
                f       = matrixFactory.f(t);
                for i = 1:Nt
                    K{i} = k1 + k2{i}; 
                end
                K = blkdiag(K{:});
                
                CD  = cell(Nt,Nt);
                PCD = cell(Nt,Nt);
                [PCD{:}] = deal(sparse(size(c,1),size(c,2)));
                
                D = exp(-1i*2*pi*((0:(Nt-1)).'*(0:(Nt-1)))/Nt);
                P = D' / Nt;
                
                for i = 1:Nt
                    for j = 1:Nt
                        CD{i,j} = omega(i)*c*D(i,j);
                    end
                end
                
                for i = 1:Nt
                    for j = 1:Nt
                        for l = 1:Nt
                            PCD{i,j} = PCD{i,j} + P(i,l) * CD{l,j};
                        end
                    end
                end
                PCD = real(cell2mat(PCD));

                x = f * 0;
                K = K + PCD;
                
                [Jm,g] = matrixFactory.G(t,x);
                r      = this.calculateResidual(x,k1,k2,c,Omega,g,f);
                x      = x(:);
                r      = r(:);
                %r      = K * x + g(:) - f;
                
                tic
                relres = norm(r) / norm(f(:))
                itter  = 0;
                while ((relres > sqrt(eps)) && (itter < 10)) || (itter == 0)
                    J  = K + blkdiag(Jm{:});
                    dy = J \ r;
                    x  = x - dy;
                    
                    x      = reshape(x,[],Nt);
                    [Jm,g] = matrixFactory.G(t,x);
                    r      = this.calculateResidual(x,k1,k2,c,Omega,g,f);
                    x      = x(:);
                    r      = r(:);
                    
                    relres = norm(r) / norm(f(:))
                    itter  = itter + 1
                end
                toc
                x = reshape(x,[],Nt);
            else
                error('');
            end
            
            x_t = fft(x,[],2);
            x_t = x_t * Omega;
            x_t = ifft(x_t,[],2);
            
            x          = mat2cell(x,size(x,1),ones(1,size(x,2)));
            x{end+1}   = x{1};
            x_t        = mat2cell(x_t,size(x_t,1),ones(1,size(x_t,2)));
            x_t{end+1} = x_t{1};
            
            %% Save the solution
            [this.X, this.X_t] = matrixFactory.doPostProcessing(x,x_t);
            solution           = Solution(this);
        end
    end
    
    methods (Static)
        %% DD Methods Here
        function r = calculateResidual(x,K1,K2,C,Omega,g,f)
            x_t = fft(x,[],2);
            x_t = x_t * Omega;
            x_t = ifft(x_t,[],2);
            r   = K1 * x + C * x_t + g - f;
            for i = 1:numel(K2)
                r(:,i) = r(:,i) + K2{i} * x(:,i);
            end
        end
        
        function [E1, E2] = extensionMatrices(I1,I2,Ib)
            E1 = [sparse(sum(I1),sum(Ib));
                  speye(sum(Ib),sum(Ib))];
            
            E2 = [sparse(sum(I2),sum(Ib));
                  speye(sum(Ib),sum(Ib))];
        end
        
        function [Mn,Mii,Mib,Mbi,Mbb] = splitMatrix(M,Ii,Ib)
            if iscell(M)
                n   = numel(M);
                Mn  = cell(n,1);
                Mii = cell(n,1);
                Mib = cell(n,1);
                Mbi = cell(n,1);
                Mbb = cell(n,1);
                for i = 1:n
                    Mii{i} = M{i}(Ii,Ii);
                    Mib{i} = M{i}(Ii,Ib);
                    Mbi{i} = M{i}(Ib,Ii);
                    Mbb{i} = M{i}(Ib,Ib);
                    Mn{i}  = [Mii{i} Mib{i};
                              Mbi{i} Mbb{i}];
                end
            else
                Mii = M(Ii,Ii);
                Mib = M(Ii,Ib);
                Mbi = M(Ib,Ii);
                Mbb = M(Ib,Ib);
                Mn  = [Mii Mib;
                       Mbi Mbb];
            end
        end
        
        function M2 = addMatrices(M1,M2)
            if iscell(M1)
                for i = 1:numel(M2)
                    M2{i} = M1{i} + M2{i};
                end
            else
                for i = 1:numel(M2)
                    M2{i} = M1 + M2{i};
                end
            end
        end
        
        function Ma = averageMatrices(M)
            N  = numel(M);
            Ma = M{1};
            for i = 2:N
                Ma = Ma + M{i};
            end
            Ma = Ma / N;
        end
        
        function [y1,y2,yb] = lowerTriangularSolve(r,I1,I2,Ib,K1ii,K1bi,C1ii,C1bi,omega,K2ii,K2bi)
            Nh = length(omega);
            
            %% Linear Region
            r1 = r(I1,:);
            r1 = fft(r1,[],2);
            y1 = 0 * r1;
            yb = zeros(sum(Ib),Nh);
            for i = 1:Nh
                y1(:,i) = (K1ii + omega(i)*C1ii) \ r1(:,i);
                yb(:,i) = K1bi * y1(:,i) + omega(i) * (C1bi * y1(:,i));
            end
            y1 = ifft(y1,[],2);
            yb = ifft(yb,[],2);
            
            %% Nonlinear Region
            r2 = r(I2,:);
            y2 = 0 * r2;
            for i = 1:Nh
                y2(:,i) = K2ii{i} \ r2(:,i);
                yb(:,i) = yb(:,i) + K2bi{i} * y2(:,i);
            end
            
            %% Interface
            yb = r(Ib,:) - yb;
        end
        
        function r = lowerTriangularMVP(y1,y2,yb,I1,I2,Ib,K1ii,K1bi,C1ii,C1bi,omega,K2ii,K2bi)
            Nh = length(omega);
            r  = zeros(length(I1),Nh);
            
            y1 = fft(y1,[],2);
            for i = 1:Nh
                r(I1,i) = K1ii * y1(:,i) + omega(i) * (C1ii * y1(:,i));
                r(Ib,i) = K1bi * y1(:,i) + omega(i) * (C1bi * y1(:,i));
            end
            r = ifft(r,[],2);
            
            for i = 1:Nh
                r(I2,i) = K2ii{i} * y2(:,i);
                r(Ib,i) = K2bi{i} * y2(:,i) + r(Ib,i);
            end
            
            r(Ib,:) = yb + r(Ib,:);
        end
        
        function y = MVP(x,K1ii,K1ib,K1bi,K1bb,C1ii,C1ib,C1bi,C1bb,omega,K2ii,K2ib,K2bi,K2bb)
            Nh = length(omega);
            x  = reshape(x,[],Nh);
            
            y  = fft(x,[],2);
            for i = 1:Nh
                z = K1ib * y(:,i) + omega(i) * (C1ib * y(:,i));
                
                z = K1ii{i,5} \ z;
                z = K1ii{i,3} * z;
                z = K1ii{i,1} \ z;
                z = K1ii{i,2} \ z;
                z = K1ii{i,4} * z;
                
                z = K1bi * z      + omega(i) * (C1bi * z);

                y(:,i) = K1bb * y(:,i) + omega(i) * (C1bb * y(:,i)) - z;
            end
            y = ifft(y,[],2);
            
            for i = 1:Nh
                z = K2ib{i} * x(:,i);
                
                z = K2ii{i,4}  * z;
                z = K2ii{i,3}' * z;
                z = K2ii{i,1}  \ z;
                z = K2ii{i,2}  \ z;
                z = K2ii{i,1}' \ z;
                z = K2ii{i,3}  * z;
                z = K2ii{i,4}  * z;
        
                z = K2bi{i} * z;
                
                y(:,i) = y(:,i) + K2bb{i} * x(:,i) - z;
            end
            
            y = y(:);
        end
        
        function y = PC(x,M1n,E1,nullSpace1)
            Nh = size(M1n,1);
            x  = reshape(x,[],Nh);
            
            y  = x;
            y  = fft(y,[],2);
            for i = 1:Nh
                z = E1 * y(:,i);
                
                z = M1n{i,5} \ z;
                z = M1n{i,3} * z;
                z = M1n{i,1} \ z;
                
                z(nullSpace1{i}) = 0;
                z = M1n{i,2} \ z;
                
                z = M1n{i,4} * z;
                
                y(:,i) = E1.' * z;
            end
            y = ifft(y,[],2);
            
            y = y(:);
        end
        
        function y = PCk(x,E,M)
            Nh = size(M,1);
            x  = reshape(x,[],Nh);
            
            y  = 0 * x;
            for i = 1:Nh
                z = E * x(:,i);
                
                z = M{i,5} \ z;
                z = M{i,3} * z;
                z = M{i,1} \ z;
                z = M{i,2} \ z;
                z = M{i,4} * z;
        
                z = E' * z;
                
                y(:,i) = z;
            end
            
            y = y(:);
        end
        
        function y = PCa(x,E,M)
            Nh = size(M,1);
            y  = reshape(x,[],Nh);
            y  = fft(y,[],2);
            
            for i = 1:Nh
                z = E * y(:,i);
                
                z = M{i,5} \ z;
                z = M{i,3} * z;
                z = M{i,1} \ z;
                z = M{i,2} \ z;
                z = M{i,4} * z;
                
                y(:,i) = E.' * z;
            end
            y = ifft(y,[],2);
            y = y(:);
        end
        
        function M = factorOverlapMatrices(K,G,C,P,omega)
            n  = length(K);
            Nh = length(omega);
            M  = cell(1,n);
            [M{:}] = deal(cell(1,5));
            for i = 1:n
                K{i} = K{i} + P{i}*G*P{i}';
                for j = 1:Nh
                    [M{i}{j,1},M{i}{j,2},M{i}{j,3},M{i}{j,4},M{i}{j,5}] = lu(K{i}+omega(j)*C{i});
                end
            end
        end
        
        function y = PCo(x,A,M,omega,Mo,E,E2T,T,R,T2E)
%             K{1} = K{1}+P{1}*G*P{1}';
%             K{2} = K{2}+P{2}*G*P{2}';
            
            Nh  = length(omega); 
            y   = 0 * x;
            
            res = reshape(x,[],Nh);
            if ~isempty(M)
                y   = y + M(res(:));
                res = x - A(y);
                res = reshape(res,[],Nh);
            end
            
            %% Reference Frame Transformation (1)
            x1 = E2T{1} * res;
            for i = 1:Nh
                x1(:,i) = T{1}{i} * x1(:,i);
            end
            x1 = T2E{1} * x1;
            
            %% Construct RHS (1)
            x1 = x1 + E{1} * res;
            
            %% Solve Frequency Domain Equations (2)
            x1 = fft(x1,[],2);
            y1 = 0 * x1;
            for i = 1:Nh
                z = x1(:,i);
                
                z = Mo{1}{i,5} \ z;
                z = Mo{1}{i,3} * z;
                z = Mo{1}{i,1} \ z;
                z = Mo{1}{i,2} \ z;
                z = Mo{1}{i,4} * z;
                
                y1(:,i) = z;
                
%                 y1(:,i) = (K{1}+omega(i)*C{1}) \ x1(:,i);
            end
            y1 = ifft(y1,[],2);
            
            %% Extend/Restrict to Boundary Unknowns (1)
            y12 = T2E{1}.' * y1;
            y1  = E{1}.'   * y1;
            for i = 1:Nh
                y12(:,i) = R{1}{i} * y12(:,i);
            end
            y12 = E2T{1}.' * y12;
            
            y1  = y1 + y12;
            
            %% Update Output
            y   = y + y1(:);
            
            %% Update Residual
            res = x - A(y);
            res = reshape(res,[],Nh);
            
            %% Reference Frame Transformation (2)
            x2 = E2T{2} * res;
            for i = 1:Nh
                x2(:,i) = T{2}{i} * x2(:,i);
            end
            x2 = T2E{2} * x2;
            
            %% Construct RHS (2)
            x2 = x2 + E{2} * res;
            
            %% Solve Frequency Domain Equations (2)
            x2 = fft(x2,[],2);
            y2 = 0 * x2;
            for i = 1:Nh
                z = x2(:,i);
                
                z = Mo{2}{i,5} \ z;
                z = Mo{2}{i,3} * z;
                z = Mo{2}{i,1} \ z;
                z = Mo{2}{i,2} \ z;
                z = Mo{2}{i,4} * z;
                
                y2(:,i) = z;
%                 y2(:,i) = (K{2} + omega(i)*C{2}) \ x2(:,i);
            end
            y2 = ifft(y2,[],2);
            
            %% Extend/Restrict to Boundary Unknowns (2)
            y22 = T2E{2}.' * y2;
            y2  = E{2}.'   * y2;
            for i = 1:Nh
                y22(:,i) = R{2}{i} * y22(:,i);
            end
            y22 = E2T{2}.' * y22;
            
            y2  = y2 + y22;
            
            %% Update Output
            y   = y  + y2(:);
        end
        
        function dx = upperTriangularSolve(y1,y2,dy,I1,I2,Ib,K1ii,K1ib,C1ii,C1ib,omega,K2ii,K2ib)
            Nh       = length(omega);
            dy       = reshape(dy,[],Nh);
            dx       = zeros(length(I1),Nh);
            dx(Ib,:) = dy;
            
            for i = 1:Nh
                z = K2ib{i} * dy(:,i);
                z = K2ii{i} \ z;
                
                y2(:,i) = y2(:,i) - z;
            end
            dx(I2,:) = y2;
            
            dy = fft(dy,[],2);
            y1 = fft(y1,[],2);
            for i = 1:Nh
                z = K1ib * dy(:,i) + omega(i) * (C1ib * dy(:,i));
                z = (K1ii + C1ii * omega(i)) \ z;
                
                y1(:,i) = y1(:,i) - z;
            end
            y1 = ifft(y1,[],2);
            dx(I1,:) = y1;
        end
        
        function varargout = factorFrequencyDomainMatrices(Ks,Cs,Kn,Cn,omega)
            Nh = length(omega);
            Ms = cell(Nh,5);
            Mn = cell(Nh,5);
            nullSpace = cell(Nh,1);
            for i = 1:Nh
                if ~isempty(Ks)
                    [Ms{i,1},Ms{i,2},Ms{i,3},Ms{i,4},Ms{i,5}] = lu(Ks+omega(i)*Cs);
                end
                if nargout > 1
                    [Mn{i,1},Mn{i,2},Mn{i,3},Mn{i,4},Mn{i,5}] = lu(Kn+omega(i)*Cn);

                    I            = (abs(diag(Mn{i,2})) < max(abs(diag(Mn{i,2}))) * sqrt(eps));
                    nullSpace{i} = I;
                    Mn{i,2}(I,:) = 0;
                    Mn{i,2}(I,I) = speye(sum(I),sum(I));
                end
            end
            varargout{1} = Ms;
            if nargout > 1
                varargout{2} = Mn;
                varargout{3} = nullSpace;
            end
        end
        
        function varargout = factorTimeDomainMatrices(Ks,Kn)
            Nh = length(Kn);
            Ms = cell(Nh,4);
            Mn = cell(Nh,5);
            for i = 1:Nh
                if ~isempty(Ks)
                    [Ms{i,1},Ms{i,2},Ms{i,3},Ms{i,4}] = ldl(Ks{i});
                end
                if nargout > 1
                    [Mn{i,1},Mn{i,2},Mn{i,3},Mn{i,4},Mn{i,5}] = lu(Kn{i});
                end
            end
            varargout{1} = Ms;
            if nargout > 1
                varargout{2} = Mn;
            end
        end
        
        function [W1,W2] = createWeightingMatrices(M1,M2)
            N  = length(M2);
            M  = length(M1);
            W1 = cell(N,1);
            W2 = cell(N,1);
            for i = 1:N
                N1 = sparse(1:M,1:M,diag(M1));
                N2 = sparse(1:M,1:M,diag(M2{i}));

%                 N1 = M1;
%                 N2 = M2{i};
                
                W1{i} = N1 + N2;
                W2{i} = W1{i}\N2;
                W1{i} = W1{i}\N1;
            end
        end
        
        %% Configure
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = HarmonicBalanceDomainDecomposition(varargin{:});
            end
        end
    end
end