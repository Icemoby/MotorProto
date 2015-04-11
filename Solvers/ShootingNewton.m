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
            %t          = linspace(0,t(end),13);
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
                Wv = blkdiag(getMatrix.PostProcessing.X2V);
                Wv = Wv * blkdiag(getMatrix.PostProcessing.Reduced2Full);
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
                
                
                %% Calculate p^th derivative estimates
                ec = zeros(1,Nt);
                ep = zeros(size(W,1),Nt);
                Cn = 0;
                for k = 1:Nt
                    h_k = t(k+1) - t(k);
                    %                         ep(:,k) = h_k * bth(1) * W * y_t{end,k};
                    %                         for i = 1:Ns
                    %                             ep(:,k) = ep(:,k) + h_k * bth(i+1) * W * y_t{i,k+1};
                    %                         end
                    %                         ep(:,k) = ep(:,k) * h_k^(-pev);
                    %ep(:,k) = (W * y_t{end,k+1} - W * y_t{end,k}) / (2 * h_k);
                    %ep(:,k) = W * y_t{end,k+1};
                    ep(:,k) = W * (y_t{end,k+1} - y_t{end,k}) / (factorial(2) * h_k);
                    ec(k)   = norm(ep(:,k));
                    Cn      = max(Cn, norm(W * y{end,k+1}));
                end
                ec = (ec ./ Cn);
                h  = diff(t);
                
                figure;
                plot(t(2:end),ec.*(h.^(pev)))
                pause(1);
                
                tol_max = (1e-3)^(pev/3);
                while max(ec.*(h.^pev)) > tol_max
                    %% Calculate New Tolerance and Time-Step Function
                    tol = mean(ec.*h.^(pev));
                    tol = max(tol, tol_max) * 2^(-pev);
                    
                 	hc = (tol./ec).^(1/pev);
                    I  = 1:(Nt/6);
                    for i = 1:5
                        hc(I) = min(hc(I), hc(I+i*Nt/6));
                    end

                    for i = 1:5
                        hc(I+i*Nt/6) = hc(I);
                    end
                    hc = [hc(end),hc];
                    
                    figure;
                    plot(reshape([t(1:(end-1));t(2:end)],1,[]), reshape([hc(2:end);hc(2:end)],1,[]));
                    
                    pause(1);
                    
                    %% Calculate New Time Points
                    hc = [hc,hc,hc];
                    t  = [t-T,t,t+T];
                	t1 = t(1:(end-1));
                    t2 = t(2:end);
                    h1 = hc(1:(end-1));
                    h2 = hc(2:end);
                    
                    s = [];
                    minax = [];
                    imin = [];
                    for i = Nt:(2*Nt)
                        j = i-1;
                        k = i+1;
                        if ((hc(i) < hc(j)*(1-sqrt(eps))) && (hc(i) < hc(k)*(1-sqrt(eps))))
                            s = cat(2,s,(t(i)+t(j))/2);
                            minax = cat(2,minax,-1);
                            imin = cat(2,imin,i);
                        end
                    end
                    I = (s > 0);
                    s = s(I);
                    I = (s <= s(1) + T/6*(1+sqrt(eps)));
                    s = s(I);
                    m = length(s)-1;
                    
                    for k = 1:m
                        if hc(imin(k)) < hc(imin(k+1)) %minax(k) == -1 %Forward from k
                            sf = [];
                            sk = s(k);
                            while sk < s(k+1)                          	
                                j = find((t1 <= sk) & (sk < t2)) + 1;
                                i = j;
                                while t(i)-sk  < hc(i)
                                    i = i + 1;
                                end
                                if t(i-1)-sk > hc(i)
                                    i = i - 1;
                                end
                                hf = min(hc(j:i));
                                sk = sk + hf;
                                sf = cat(2,sf,sk);
                            end
                            sf = (sf-s(k)) * (s(k+1)-s(k)) / (sf(end)-s(k)) + s(k);
                            sf(end) = [];
                            s = cat(2,s,sf);
                        else %backward from k+1
                            sb = [];
                            sk = s(k+1);
                            while sk > s(k)                            
                                j = find((t1 < sk) & (sk <= t2)) + 1;
                                i = j;
                                while sk-t(i-1) < hc(i)
                                    i = i - 1;
                                end
                                if sk-t(i) > hc(i)
                                    i = i + 1;
                                end
                                hb = min(hc(i:j));
                                sk = sk - hb;
                                sb = cat(2,sb,sk);
                            end
                            sb = (sb - s(k+1)) * (s(k+1) - s(k)) / (s(k+1)-sb(end)) + s(k+1);
                            sb(end) = [];
                            s = cat(2,s,sb);
                        end
                    end
                    s = sort(s);
                    s(end) = [];
                    t = [s,s+T/6,s+T/3,s+T/2,s+2*T/3,s+5*T/6];
                    t = mod(t,T);
                    t = sort(t);
                    t = [t,t(1)+T];
                    %%
%                     I6 = find((t>0) & (t<=T/6));
%                     [hmin,imin] = min(hc(I6));
%                     tmin        = t(I6(imin));
%                     
%                     t1 = t(1:(end-1));
%                     t2 = t(2:end);
%                     h1 = hc(1:(end-1));
%                     h2 = hc(2:end);
%                     
%                     sb = tmin;
%                     sk = tmin;
% %                     hk = hmin;
%                     while sb(end) > 0
%                         j = find((t1 < sk) & (sk <= t2)) + 1;
%                         i = j;
%                         while sk-t(i-1) < hc(i)
%                             i = i - 1;
%                         end
%                         
%                         if sk-t(i) > hc(i)
%                             i = i + 1;
%                         end
%                         
%                         hmin = min(hc(i:j));
%                         
%                         sk   = sb(end) - hmin;
%                         sb   = cat(2,sb,sk);
%                         
%                         stest = sb(end) + (t1.*h2-t2.*h1) ./ (t2-t1);
%                         stest = stest ./ (1+(h2-h1)./(t2-t1));
%                         valid = (stest >= t1) & (stest <= t2);
%                         stest = stest(valid);
%                         stest = stest(stest < sb(end));
%                         [~,i] = min(abs(stest-sb(end)));
%                         stest = stest(i);
%                         htest = sb(end)-stest;
%                         
%                         i = find(t>=stest,   1, 'first');
%                         j = find(t<=sb(end), 1, 'last');
%                         
%                         if i > j
%                             hmin = inf;
%                         else
%                         	hmin = min(hc(i:j));
%                         end
%                         hmin = min(hmin,htest);
%                         hmin = min(hmin,hk);
%                         
%                         sk = sb(end) - hmin;
%                         sb = cat(2,sb,sk);
%                         
%                         i = find(t<=sk, 1, 'last');
%                         j = find(t>=sk, 1, 'first');
%                         hk = (hc(j)*(sk-t(i))+hc(i)*(t(j)-sk))/(t(j)-t(i));
%                     end
%                     sb = fliplr(sb);
%                     
%                     if length(sb) > 1
%                         sb = (sb-sb(1)) / (sb(end)-sb(1)) * tmin;
%                     end
%                     sb(1) = 0;
% 
%                     I6 = find((t>0) & (t<=T/6));
%                     [hmin,imin] = min(hc(I6));
%                     tmin        = t(I6(imin));
%                     
%                     sf = tmin;
%                     sk = tmin;
% %                     hk = hmin;
%                     while sf(end) < T/6
%                         j = find((t1 < sk) & (sk <= t2)) + 1;
%                         i = j;
%                         while t(i)-sk  < hc(i)
%                             i = i + 1;
%                         end
%                         
%                         if t(i-1)-sk > hc(i)
%                             i = i - 1;
%                         end
%                         
%                         hmin = min(hc(j:i));
%                         sk   = sf(end) + hmin;
%                         sf   = cat(2,sf,sk);
%                         stest = sf(end) + (t2.*h1-t1.*h2) ./ (t2-t1);
%                         stest = stest ./ (1-(h2-h1)./(t2-t1));
%                         valid = (stest >= t1) & (stest <= t2);
%                         stest = stest(valid);
%                         stest = stest(stest > sf(end));
%                         [~,i] = min(abs(stest-sf(end)));
%                         stest = stest(i);
%                         htest = stest-sf(end);
%                         
%                         i = find(t>=sf(end), 1, 'first');
%                         j = find(t<=stest,   1, 'last');
%                         
%                         if i > j
%                             hmin = inf;
%                         else
%                             hmin = min(hc(i:j));
%                         end
%                         hmin = min(hmin,htest);
%                         hmin = min(hmin,hk);
%                         
%                         sk = sf(end) + hmin;
%                         sf = cat(2,sf,sk);
%                         
%                         i = find(t<=sk, 1, 'last');
%                         j = find(t>=sk, 1, 'first');
%                         hk = (hc(j)*(sk-t(i))+hc(i)*(t(j)-sk))/(t(j)-t(i));
%                     end
%                     if length(sf) > 1
%                         sf = (sf-tmin) / (sf(end)-tmin) * (T/6-tmin) + tmin;
%                     end
%                     sf(1) = [];
%                     
%                     s = sort([sb,sf]);
%                     if abs(s(end)-T/6) < sqrt(eps)*T/6
%                         s(end) = [];
%                     end
%                     
%                     t = [s,s+T/6,s+T/3,s+T/2,s+T*2/3,s+T*5/6];
%                     t = mod(t,T);
%                     t = sort(t);
%                     t = [t,t(1)+T];
                    %%
                    hold on;
                    plot(reshape([t(1:end-1);t(2:end)],1,[]), reshape([diff(t);diff(t)],1,[]),':x');
                    pause(1);
                    
                    length(t)
                    
                    %% Initialize
                    Nt = length(t) - 1;
                    y  = repmat(y(:,1),1,Nt+1);
                    y_t = repmat(y_t(:,1),1,Nt+1);
                    this.Times = t;
                    this.TimePoints = Nt;
                    
                    %% Solve
                    if tol > tol_max
                        this.MaxShootingIterations = 10;
                        this.ShootingTolerance = 1e-2;
                    else
                        this.MaxShootingIterations = 100;
                        this.ShootingTolerance = 1e-6;
                    end
                    [y,y_t,t] = this.solve_stored(y, y_t, t, c, d, p, q, Nt, Ns, getMatrix);
                
                    
                    %% Calculate p^th derivative estimates
                    ec = zeros(1,Nt);
                    ep = zeros(size(W,1),Nt);
                    Cn = 0;
                    for k = 1:Nt
                        h_k = t(k+1) - t(k);
%                         ep(:,k) = h_k * bth(1) * W * y_t{end,k};
%                         for i = 1:Ns
%                             ep(:,k) = ep(:,k) + h_k * bth(i+1) * W * y_t{i,k+1};
%                         end
%                         ep(:,k) = ep(:,k) * h_k^(-pev);
                        %ep(:,k) = (W * y_t{end,k+1} - W * y_t{end,k}) / (2 * h_k);
                        %ep(:,k) = W * y_t{end,k+1};
                        ep(:,k) = W * (y_t{end,k+1} - y_t{end,k}) / (factorial(2) * h_k);
                        ec(k)   = norm(ep(:,k));
                        Cn      = max(Cn, norm(W * y{end,k+1}));
                    end
                    ec = (ec ./ Cn);
                    h  = diff(t);
                    figure;
                    plot(t(2:end),ec.*(h.^(pev)))
                    pause(1);
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