classdef ShootingNewton < Solver
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
    % See also MotorProto, Solver
    
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
        ShootingTolerance     = 1e-6;
        NewtonTolerance       = 1e-6;
        GMRESTolerance        = 1e-6;
        MaxShootingIterations = 100;
        MaxNewtonIterations   = 100;
        MaxGMRESIterations    = 100;
        RungeKuttaStages      = 3;
        SymmetricJacobian     = false;
    end
    
    methods
        %% Constructor
        function this = ShootingNewton(varargin)
            if nargin > 0
                for iArg = 1:2:nargin
                    this.(varargin{iArg}) = varargin{iArg+1};
                end
            end
        end
        
        %% Solve
        function solution = solve(this, model, y0)
            %% Setup Matrices
            getMatrix = DynamicMatrixFactory(copy(model));
            this.Matrices = getMatrix;
            
            Nx = length(getMatrix.f(0));
            
            %% Get RK-Method
            Ns = this.RungeKuttaStages;
            [~,~,c,d,p,q,bu,bh,bth] = this.getButcherTable(Ns);
            
            %% Initial Condition
            t          = model.getTimePoints(this.TimePoints);
            %t          = linspace(0,t(end),37);
            this.Times = t;
            Nt         = numel(t) - 1;
         	y          = cell(Ns, Nt+1);
            y_t        = cell(Ns, Nt+1);
            [y{:}]     = deal(zeros(Nx, 1));
            [y_t{:}]   = deal(zeros(Nx, 1));
            
            if nargin < 3 || isempty(y0)
                y{end,1}  = zeros(Nx, 1);
            elseif iscell(y0)
                assert(length(y0{end, 1}) == Nx);
            	y{end,1} = y0{end, 1};
            else
                assert(length(y0) == Nx);
                y{end,1} = y0;
            end
            
            %% Begin Simulation            
            if this.Verbose
                display(sprintf('\nShooting-Newton %d/%d',Ns,Nt));
            end
            
            tic;
            if this.StoreDecompositions
                [y,y_t,t] = this.solve_stored(y, y_t, t, c, d, p, q, Nt, Ns, getMatrix);
                T = t(end) - t(1);
                
                %% Error Indicator
                Wa = blkdiag(getMatrix.PostProcessing.SobolevA);
             	Wa = Wa * blkdiag(getMatrix.PostProcessing.Reduced2Full);
%                 Wv = blkdiag(getMatrix.PostProcessing.X2V);
%                 Wv = Wv * blkdiag(getMatrix.PostProcessing.Reduced2Full);
%                 ev = this.rkError(t,y(:,2:end),y_t(:,2:end),bh,Wv);
%                 ea = this.rkError(t,y(:,2:end),y_t(:,2:end),bh,Wa);
                W = Wa;
%                 pv  = 2;
%                 rv  = 2;
                pev = 2;
%                 
%                 pa  = 3;
%                 ra  = 2;
%                 pea = 1;
                
                tol = inf;
                tol_max = 1e-3;
                while tol > tol_max
                    %% Calculate p^th derivative estimates
                    ec = zeros(1,Nt);
                    ep = zeros(size(W,1),Nt);
                    Cn = 0;
                    for k = 1:Nt
                        h_k = t(k+1) - t(k);
                        ep(:,k) = h_k * bth(1) * W * y_t{end,k};
                        for i = 1:Ns
                            ep(:,k) = ep(:,k) + h_k * bth(i+1) * W * y_t{i,k+1};
                        end
                        ep(:,k) = ep(:,k) * h_k^(-pev);
                        ec(k)   = norm(ep(:,k));
                        Cn      = max(Cn, norm(W * y{end,k+1}));
                    end
                    ec = (ec ./ Cn);
                    
                    %% Calculate New Tolerance and Time-Step Function
                    h   = diff(t);
                    tol = (sum(h.^2)^(pev)) * (sum(h.*ec.^(-1/pev)))^(-pev);
                    tol = tol * 2^(-pev)
                    if tol < tol_max
                        tol = tol_max;
                    end
                 	hc  = (tol./ec).^(1/pev);
                    hmin = min(hc);
                    hmax = max(hc);
                    
                    %% Do Harmonic Linear Interpolation
                    N  = 12 * Nt / 2;
                    T  = t(end);
                    t1 = t(1:(end-1));
                    t2 = t(2:end);
                    
                    w = 2*pi*[6:6:N fliplr(-(6:6:N))] / T;

                    M1 = [ (t2-t1)/2;
                            -(exp(1i*w.'*t1)-exp(1i*w.'*t2)-1i*(w.'*t1).*exp(1i*w.'*t1)+1i*(w.'*t2).*exp(1i*w.'*t1))./((w.').^2*(t1-t2))];
                    M2 = [ (t2-t1)/2;
                            (exp(1i*w.'*t1)-exp(1i*w.'*t2)-1i*(w.'*t1).*exp(1i*w.'*t2)+1i*(w.'*t2).*exp(1i*w.'*t2))./((w.').^2*(t1-t2))];
                    M  = M1 + circshift(M2,[0,1]);
                    M  = M / T;
                    w = 2*pi*[0:6:N fliplr(-(6:6:N))] / T;

                    HC = (M * hc.').';
                    
                    s = linspace(0,t(end),1000);
                    test = zeros(size(s));
                    for i = 1:length(test)
                        M = exp(-1i*w.'*s(i));
                        test(i) = real(HC * M);
                    end
                    figure;plot(s,test);
                    hold on;scatter(t(1:(end-1)),hc);
                    figure;plot(t(2:end),ec.*diff(t).^(pev))
                    pause(1);
                    
                    %% Calculate New Time Points
                    temp = toc;
                    s = t(end);
                    while s(end) > T*5/6
                        Mend = exp(-1i*w.'*s(end)) ./ (-1i*w.');
                        Mend(1) = s(end);
                        Hend = real(HC * Mend);
                        
                        sk = s(end) - hmax;
                        Mk = exp(-1i*w.'*sk) ./ (-1i*w.');
                        Mk(1) = sk;
                        Hk = real(HC * Mk);
                        
                        r = Hend-Hk-(s(end)-sk)^2;
                        while abs(r) > sqrt(eps) * hmin
                            mk = exp(-1i*w.'*sk);
                            hk = real(HC * mk);
                            
                            J  = 2*(s(end)-sk)-hk;
                            sk = sk - J\r;
                            
                            Mk = exp(-1i*w.'*sk) ./ (-1i*w.');
                            Mk(1) = sk;
                            Hk = real(HC * Mk);
                            
                            r = Hend-Hk-(s(end)-sk)^2;
                        end
                        s(end+1) = sk;
%                         mk  = exp(-1i*w.'*(s(end)-eps));
%                         hk  = real(HC * mk);
%                         dmk = -1i*w.'.*exp(-1i*w.'*(s(end)-eps));
%                         dhk = real(HC * dmk);
%                         
%                         if dhk < 0
%                             snew = s(end) - hk;
%                         else
%                             snew = s(end) - hmin;
%                             
%                             mk  = exp(-1i*w.'*snew);
%                             hk  = real(HC * mk);
%                             dmk = -1i*w.'.*exp(-1i*w.'*snew);
%                             dhk = real(HC * dmk);
%                             r = snew + hk - s(end);
%                             J = dhk + 1;
%                             
%                             while r > sqrt(eps) * hmin
%                                 snew = snew - 0.5*(J \ r);
%                                 mk  = exp(-1i*w.'*snew);
%                                 hk  = real(HC * mk);
%                                 dmk = -1i*w.'.*exp(-1i*w.'*snew);
%                                 dhk = real(HC * dmk);
%                                 r = snew + hk - s(end);
%                                 J = dhk + 1;
%                             end
%                         end
%                         s(end+1) = snew;
                    end
                    toc - temp
                    
                    s = fliplr(s);
                    s = s - T * 5 /6;
                    s = s - s(1) * (6*s - T) / (6*s(1) - T);
                    s(1) = 0;
                    s(end) = [];
                    t = [s,s+T/6,s+T/3,s+T/2,s+T*2/3,s+T*5/6,T];
                    length(t)
                    
                    %% Initialize
                    Nt = length(t) - 1;
                    y  = repmat(y(:,1),1,Nt+1);
                    y_t = repmat(y_t(:,1),1,Nt+1);
                    this.Times = t;
                    this.TimePoints = Nt;
                    
                    %% Solve
                    if tol > tol_max
                        this.MaxShootingIterations = 0;
                    else
                        this.MaxShootingIterations = 100;
                    end
                    [y,y_t,t] = this.solve_stored(y, y_t, t, c, d, p, q, Nt, Ns, getMatrix);
                end
            else
            	[y,y_t,t] = this.solve_unstored(y, y_t, t, c, d, p, q, Nt, Ns, getMatrix);
            end
            this.SimulationTime = toc;
            
           	if this.Verbose
                display(sprintf('\nSimulation Time = %0.3g seconds', this.SimulationTime));
            end
            
            %% Post Processing
            Nt  = length(t) - 1;
            x   = cell(1, Nt+1);
            x_t = cell(1, Nt+1);
            for k = 1:Nt
                x{k+1}   = y{end,k+1};
                x_t{k+1} = y_t{end,k+1};
            end
            x{1}   = x{end};
            x_t{1} = x_t{end};
            
            [x, x_t] = getMatrix.doPostProcessing(x, x_t);
            this.Y   = y;
            this.Y_t = y_t;
            this.X   = x;
            this.X_t = x_t;
            solution = Solution(this);
        end
        
        function [y,y_t,t] = solve_unstored(this, y, y_t, t, c, d, p, q, Nt, Ns, getMatrix)
            %% Algorithm Parameters
          	maxShooting = this.MaxShootingIterations;
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            
            shootingTol = this.ShootingTolerance;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            
            verbose     = this.Verbose;
            
            %% Minimal Matrix Storage
            Ca = cell(1, Ns);
            
            nShooting   = 0;
            shootingRes = 1;
            while shootingRes > shootingTol && nShooting <= maxShooting
                nShooting = nShooting + 1;
                for k = 1:Nt
                    h_k = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * h_k;
                        h_ik = h_k / d(i,i);
                        
                        y{i,k+1} = y{Ns,k};
                        
                        [~, rpq] = getMatrix.G(t(k), y{Ns,k}, h_ik);
                        Kpq = getMatrix.K(t(k), h_ik);
                        Cpq = getMatrix.C(t(k), h_k / p(i), h_ik);
                        rpq = rpq + Kpq * y{Ns,k};
                        rpq = rpq - getMatrix.f(t(k), h_ik);
                        rpq = q(i) * rpq - Cpq * y{Ns,k};
                        
                        for j = 1:i
                            Ca{j} = getMatrix.C(t_ik, h_k / d(i,j), h_ik);
                        end
                        
                        f_ik = getMatrix.f(t_ik, h_ik);
                        K_ik = getMatrix.K(t_ik, h_ik);
                        
                        %% Newton Iteration
                        nNewton   = 0;
                        newtonRes = 1;
                        while newtonRes > newtonTol && nNewton < maxNewton
                            nNewton = nNewton + 1;
                            [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k+1}, h_ik);
                            
                            J_ik = K_ik + Ca{i} + G_ik;
                            r_ik = K_ik * y{i,k+1} + g_ik - f_ik + rpq;
                            for j = 1:i
                                r_ik = r_ik + Ca{j} * y{j,k+1};
                            end
                            
                            newtonRes = norm(r_ik) / norm(g_ik - f_ik);
                            
                            r_ik     = J_ik \ r_ik;
                            y{i,k+1} = y{i,k+1} - r_ik;
                        end
                        
                        %% Calculate Stage Derivative
                        y_t{i, k+1} = (-p(i)/h_k) * y{Ns,k} - q(i)*y_t{Ns,k};
                        for j = 1:i
                            y_t{i,k+1} = y_t{i,k+1} + (d(i,j)/h_k)*y{j,k+1};
                        end
                    end
                end
                
                ei = this.rkError(t,y,y_t,bh,bth,2:Nt,getMatrix);
                figure;
                plot(ei);
                
                r = y{end,end} - y{end,1};
                shootingRes = norm(r) / norm(y{end,end});
                
                if verbose
                    display(sprintf('\nIteration %d, Residual = %0.3g', nShooting, shootingRes));
                end
                
                if (shootingRes > shootingTol) && (maxGMRES > 0)
                    A = @ShootingNewton.MVPUnstored;
                    if verbose
                        dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, y, t, c, d, p, q, Ns, Nt, getMatrix);
                    else
                        [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, y, t, c, d, p, q, Ns, Nt, getMatrix);
                    end
                    y{end,1} = y{end,1} - dy;
                elseif maxGMRES == 0
                    y{end,1} = y{end,end};
                    display(sprintf(' '));
                end
            end
        end

        function [y,y_t,t] = solve_stored(this, y, y_t, t, c, d, p, q, Nt, Ns, getMatrix)
            %% Algorithm Parameters
          	maxShooting = this.MaxShootingIterations;
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            
            shootingTol = this.ShootingTolerance;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            
            isSymmetric = this.SymmetricJacobian;
            
            verbose     = this.Verbose;

            %% Store Everything
            K   = cell(Ns, Nt);
            f   = cell(Ns, Nt);
            Ca  = cell(Ns, Ns, Nt);
            Jpq = cell(Ns, Nt);
            if this.SymmetricJacobian
                L = cell(Ns, Nt);
                D = cell(Ns, Nt);
                S = cell(Ns, Nt);
                P = cell(Ns, Nt);
            else
                L = cell(Ns, Nt);
                U = cell(Ns, Nt);
                P = cell(Ns, Nt);
                Q = cell(Ns, Nt);
                R = cell(Ns, Nt);
            end
            
            for k = 1:Nt
                h_k  = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * h_k;
                    h_ik = h_k / d(i,i);
                    for j = 1:i
                        Ca{i,j,k} = getMatrix.C(t_ik, h_k / d(i,j), h_ik);
                    end
                    K{i,k}  = getMatrix.K(t_ik, h_ik);
                    f{i,k}  = getMatrix.f(t_ik, h_ik);
                end
            end
            
            %% Start
            nShooting   = 0;
            shootingRes = 1;
            while shootingRes > shootingTol && nShooting <= maxShooting
                nShooting = nShooting + 1;
                for k = 1:Nt
                    h_k  = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * h_k;
                        h_ik = h_k / d(i,i);
                        
                    	y{i, k+1} = y{Ns, k};
                        
                        [Gpq, rpq] = getMatrix.G(t(k), y{Ns,k}, h_ik);
                        Kpq         = getMatrix.K(t(k), h_ik);
                        Cpq         = getMatrix.C(t(k), h_k / p(i), h_ik);
                        Jpq{i,k}    = Cpq - q(i) * (Gpq + Kpq);
                        rpq         = rpq + Kpq * y{Ns,k};
                        rpq         = rpq - getMatrix.f(t(k), h_ik);
                        rpq         = q(i) * rpq - Cpq * y{Ns,k};
                        
                        %% Newton Iteration
                        nNewton   = 0;
                        newtonRes = 1;
                        while newtonRes > newtonTol && nNewton < maxNewton
                            nNewton = nNewton + 1;
                            
                            [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k+1}, h_ik);
                            
                            r_ik = K{i,k} * y{i,k+1} + g_ik - f{i,k} + rpq;
                            for j = 1:i
                                r_ik = r_ik + Ca{i,j,k} * y{j,k+1};
                            end
                            
                            newtonRes = norm(r_ik) / norm(g_ik - f{i,k});
                            
                            J_ik = K{i,k} + Ca{i,i,k} + G_ik;
                            if isSymmetric
                                [L_ik, D_ik, P_ik, S_ik] = ldl(J_ik);
                                r_ik = S_ik  * r_ik;
                                r_ik = P_ik' * r_ik;
                                r_ik = L_ik  \ r_ik;
                                r_ik = D_ik  \ r_ik;
                                r_ik = L_ik' \ r_ik;
                                r_ik = P_ik  * r_ik;
                                r_ik = S_ik  * r_ik;
                            else
                                [L_ik, U_ik, P_ik, Q_ik, R_ik] = lu(J_ik);
                                r_ik = R_ik \ r_ik;
                                r_ik = P_ik * r_ik;
                                r_ik = L_ik \ r_ik;
                                r_ik = U_ik \ r_ik;
                                r_ik = Q_ik * r_ik;
                            end
                            y{i,k+1} = y{i,k+1} - r_ik;
                        end
                        
                        %% Calcualte Stage Derivative
                        y_t{i, k+1} = (-p(i)/h_k) * y{Ns,k} - q(i)*y_t{Ns,k};
                        for j = 1:i
                            y_t{i,k+1} = y_t{i,k+1} + (d(i,j)/h_k)*y{j,k+1};
                        end
                        
                        if isSymmetric
                            L{i,k} = L_ik;
                            D{i,k} = D_ik;
                            P{i,k} = P_ik;
                            S{i,k} = S_ik;
                        else
                            L{i,k} = L_ik;
                            U{i,k} = U_ik;
                            P{i,k} = P_ik;
                            Q{i,k} = Q_ik;
                            R{i,k} = R_ik;
                        end
                    end
                end
                
                r = y{end,end} - y{end,1};
                shootingRes = norm(r) / norm(y{end, end});
                
                if verbose
                    display(sprintf('\nIteration %d, Residual = %0.3g', nShooting, shootingRes));
                end

                if (shootingRes > shootingTol) && (maxGMRES > 0)
                    if isSymmetric
                        A = @ShootingNewton.MVPStoredLDL;
                        if verbose
                            dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Jpq, L, D, P, S, Ns, Nt);
                        else
                            [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, K, Ca, Jpq, L, D, P, S, Ns, Nt);
                        end
                    else
                        A = @ShootingNewton.MVPStoredLU;
                        if verbose
                            dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Jpq, L, U, P, Q, R, Ns, Nt);
                        else
                            [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Jpq, L, U, P, Q, R, Ns, Nt);
                        end
                    end
                    y{end,1}   = y{end,1} - dy;
                    y_t{end,1} = y_t{end,end};
                elseif maxGMRES == 0
                    y{end,1}   = y{end,end};
                    y_t{end,1} = y_t{end,end};
                    display(sprintf(' '));
                end
            end
        end
    end
    
    methods (Static)
        %% Stored Matrices
        function mvp = MVPStoredLDL(z_s0, Ca, Jpq, L, D, P, S, Ns, Nt)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                for i = 1:Ns         
                    rhs = Jpq{i,k} * z_sk;
                    for j = 1:(i - 1);
                        rhs = rhs - Ca{i,j,k} * z{j};
                    end
                    z{i} = S{i,k} * (P{i,k} * (L{i,k}' \ (D{i,k} \ (L{i,k} \ (P{i,k}' * (S{i,k} * rhs))))));
                end
                z_sk = z{Ns};
            end
            mvp = z_sk - z_s0;
        end
        
        function mvp = MVPStoredLU(z_s0, Ca, Jpq, L, U, P, Q, R, Ns, Nt)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                for i = 1:Ns         
                    r_ik = Jpq{i,k} * z_sk;
                    for j = 1:(i - 1);
                        r_ik = r_ik - Ca{i,j,k} * z{j};
                    end
                    z{i} = Q{i,k} * (U{i,k} \ (L{i,k} \ (P{i,k} * (R{i,k} \ r_ik))));
                end
                z_sk = z{Ns};
            end
            mvp = z_sk - z_s0;
        end
        
        %% Unstored Matrices
        function mvp = MVPUnstored(z_s0, y, t, c, d, p, q, Ns, Nt, getMatrix)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                h_k = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * h_k;
                    h_ik = h_k / d(i,i);
                    
                    Jpq = getMatrix.G(t(k), y{Ns,k}, h_ik);
                    Jpq = Jpq + getMatrix.K(t(k), h_ik);
                    Jpq = getMatrix.C(t(k), h_k / p(i), h_ik) - q(i)*Jpq;
                    r_ik = Jpq * z_sk;
                    for j = 1:(i-1)
                        C    = getMatrix.C(t_ik, h_k / d(i,j), h_ik);
                        r_ik = r_ik - C * z{j};
                    end

                  	J_ik = getMatrix.K(t_ik, h_ik);
                    J_ik = J_ik + getMatrix.G(t_ik, y{i,k+1}, h_ik);
                    J_ik = J_ik + getMatrix.C(t_ik, h_ik, h_ik);
                    
                    z{i} = J_ik \ r_ik;
                end
                z_sk = z{Ns};                
            end
            mvp = z_sk - z_s0;
        end
      
        %% Configure
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = ShootingNewton(varargin{:});
            end
        end
    end
end