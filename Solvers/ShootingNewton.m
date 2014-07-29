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
        ShootingTolerance     = 1e-6;
        NewtonTolerance       = 1e-6;
        GMRESTolerance        = 1e-6;
        MaxShootingIterations = 100;
        MaxNewtonIterations   = 100;
        MaxGMRESIterations    = 100;
        RungeKuttaStages      = 2;
        StorageLevel          = 0;
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
            Ns        = this.RungeKuttaStages;
            [A, ~, c] = this.getButcherTable(Ns);
            alpha     = inv(A);
            gamma     = sum(alpha,2);
            
            %% Initial Condition
            t          = model.getTimePoints(this.TimePoints);
            this.Times = t;
            Nt         = numel(t) - 1;
         	y          = cell(Ns, Nt+1);
            [y{:}]     = deal(zeros(Nx, 1));
            
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
                display(sprintf('GMRES %d/%d\n',Ns,Nt));
            end
            
            tic;
            switch this.StorageLevel
                case 0
                    y = this.solve_level_0(y, t, alpha, gamma, c, Nt, Ns, getMatrix);
                case 1
                    error('#TODO');
                case 2
                    y = this.solve_level_2(y, t, alpha, gamma, c, Nt, Ns, getMatrix);
                case 3
                    y = this.solve_level_3(y, t, alpha, gamma, c, Nt, Ns, getMatrix);
                otherwise
                    error('Storage level must be between 0 and 3');
            end
            this.SimulationTime = toc;
            
           	if this.Verbose
                display(sprintf('Simulation Time = %0.3g seconds\n', this.SimulationTime));
            end
            
            %% Post Processing
            x   = cell(1, Nt+1);
            x_t = cell(1, Nt+1);
            for k = 1:Nt
                hk = t(k+1) - t(k);
                
                x{k} = y{end,k};
                
                x_t{k+1} = (-gamma(end) / hk) * y{end,k};
                for i = 1:Ns
                    x_t{k+1} = x_t{k+1} + (alpha(end,i) / hk) * y{i,k+1};
                end
            end
            x{end} = x{1};
            x_t{1} = x_t{end};
            y(:,1) = y(:,end);
            
            [x, x_t] = getMatrix.doPostProcessing(x, x_t);
            this.Y   = y;
            this.X   = x;
            this.X_t = x_t;
            solution = Solution(this);
        end
        
        function y = solve_level_0(this, y, t, alpha, gamma, c, Nt, Ns, getMatrix)
            %% Algorithm Parameters
          	maxShooting = this.MaxShootingIterations;
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            
            shootingTol = this.ShootingTolerance;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            
            verbose     = this.Verbose;
            
            %% Allocate Matrix Storage
            Ca = cell(1, Ns);
            d_ik = y{end, 1} * 0;

            nShooting   = 0;
            shootingRes = 1;
            while shootingRes > shootingTol && nShooting <= maxShooting
                nShooting = nShooting + 1;
                for k = 1:Nt
                    hk  = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        if i == 1
                            y{1, k+1} = y{Ns, k} + hk * c(i) * d_ik;
                        else
                            y{i, k+1} = y{i-1, k+1} + hk * (c(i) - c(i-1)) * d_ik;
                        end
                        
                        for j = 1:i
                            Ca{j} = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                        end
                        Cg   = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                        f_ik = getMatrix.f(t_ik, h_ik);
                        K_ik = getMatrix.K(t_ik, h_ik);
                        
                        %% Perform newton - raphson itteration
                        nNewton   = 0;
                        newtonRes = 1;
                        while newtonRes > newtonTol && nNewton < maxNewton
                            nNewton = nNewton + 1;
                            [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k+1}, h_ik);
                            
                            J_ik = K_ik + Ca{i} + G_ik;
                            
                            r_ik = K_ik * y{i,k+1} + g_ik - f_ik;
                            for j = 1:i
                                r_ik = r_ik + Ca{j} * y{j,k+1};
                            end
                            r_ik = r_ik - Cg * y{Ns, k};
                            
                            newtonRes = norm(r_ik) / norm(g_ik - f_ik);
                            
                            r_ik = J_ik \ r_ik;
                            y{i,k+1} = y{i,k+1} - r_ik;
                        end
                        
                        d_ik = (-gamma(i) / hk) * y{end,k};
                        for j = 1:i
                            d_ik = d_ik + (alpha(i,j) / hk) * y{j,k+1};
                        end
                    end
                end
                
                r = y{end,end} - y{end,1};
                shootingRes = norm(r) / norm(y{end,end});
                
                if verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', nShooting, shootingRes));
                end
                
                if shootingRes > shootingTol && maxGMRES > 0
                    A = @ShootingNewton.MVP0;
                    if verbose
                        dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, y, t, alpha, gamma, c, Ns, Nt, getMatrix);
                    else
                        [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, y, t, alpha, gamma, c, Ns, Nt, getMatrix);
                    end
                    y{end,1} = y{end,1} - dy;
                else
                    y{end,1} = y{end,end};
                    display(sprintf(' '));
                end
            end
        end
        
        function y = solve_level_2(this, y, t, alpha, gamma, c, Nt, Ns, getMatrix)
            %% Algorithm Parameters
          	maxShooting = this.MaxShootingIterations;
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            
            shootingTol = this.ShootingTolerance;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            
            verbose     = this.Verbose;
            
            %% Allocation Level 3, Store Everything
            K  = cell(Ns, Nt);
            f  = cell(Ns, Nt);
            Ca = cell(Ns, Ns, Nt);
            Cg = cell(Ns, Nt);
            J  = cell(Ns, Nt);
            
            for k = 1:Nt
                hk  = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    for j = 1:i
                        Ca{i,j,k} = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                    end
                    Cg{i,k} = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                    K{i,k}  = getMatrix.K(t_ik, h_ik);
                    f{i,k}  = getMatrix.f(t_ik, h_ik);
                end
            end
            d_ik = f{1,1} * 0;
            
            %% Start
            nShooting   = 0;
            shootingRes = 1;
            while shootingRes > shootingTol && nShooting <= maxShooting
                nShooting = nShooting + 1;
                for k = 1:Nt
                    hk  = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i, i);
                        
                        if i == 1
                            y{1, k+1} = y{Ns, k} + hk * c(i) * d_ik;
                        else
                            y{i, k+1} = y{i-1, k+1} + hk * (c(i) - c(i-1)) * d_ik;
                        end
                        
                        nNewton   = 0;
                        newtonRes = 1;
                        while newtonRes > newtonTol && nNewton < maxNewton
                            nNewton = nNewton + 1;
                            
                            [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k+1}, h_ik);
                            
                            r_ik = K{i,k} * y{i,k+1} + g_ik - f{i,k};
                            for j = 1:i
                                r_ik = r_ik + Ca{i,j,k} * y{j,k+1};
                            end
                            r_ik = r_ik - Cg{i,k} * y{Ns, k};
                            
                            newtonRes = norm(r_ik) / norm(g_ik - f{i,k});
                            
                            J{i,k}   = K{i,k} + Ca{i,i,k} + G_ik;
                            r_ik     = J{i,k} \ r_ik;
                            
                            y{i,k+1} = y{i,k+1} - r_ik;
                        end
                        
                        d_ik = (-gamma(i) / hk) * y{end,k};
                        for j = 1:i
                            d_ik = d_ik + (alpha(i,j) / hk) * y{j,k+1};
                        end
                    end
                end
                
                r = y{end,end} - y{end,1};
                shootingRes = norm(r) / norm(y{end, end});
                
                if verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', nShooting, shootingRes));
                end

                if shootingRes > shootingTol
                    A = @ShootingNewton.MVP2;
                    if verbose
                        dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Cg, J, Ns, Nt);
                        display(sprintf(' '));
                    else
                        [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Cg, J, Ns, Nt);
                    end
                    y{end,1} = y{end,1} - dy;
                end
            end
        end
        
        function y = solve_level_3(this, y, t, alpha, gamma, c, Nt, Ns, getMatrix)
            %% Algorithm Parameters
          	maxShooting = this.MaxShootingIterations;
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            
            shootingTol = this.ShootingTolerance;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            
            isSymmetric = this.SymmetricJacobian;
            
            verbose     = this.Verbose;
            
            %% Allocation Level 3, Store Everything
            K  = cell(Ns, Nt);
            f  = cell(Ns, Nt);
            Ca = cell(Ns, Ns, Nt);
            Cg = cell(Ns, Nt);
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
                hk  = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    for j = 1:i
                        Ca{i,j,k} = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                    end
                    Cg{i,k} = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                    K{i,k}  = getMatrix.K(t_ik, h_ik);
                    f{i,k}  = getMatrix.f(t_ik, h_ik);
                end
            end
            d_ik = f{1,1} * 0;
            
            %% Start
            nShooting   = 0;
            shootingRes = 1;
            while shootingRes > shootingTol && nShooting <= maxShooting
                nShooting = nShooting + 1;
                for k = 1:Nt
                    hk  = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i, i);
                        
                        if i == 1
                            y{1, k+1} = y{Ns, k} + hk * c(i) * d_ik;
                        else
                            y{i, k+1} = y{i-1, k+1} + hk * (c(i) - c(i-1)) * d_ik;
                        end
                        
                        nNewton   = 0;
                        newtonRes = 1;
                        while newtonRes > newtonTol && nNewton < maxNewton
                            nNewton = nNewton + 1;
                            
                            [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k+1}, h_ik);
                            
                            r_ik = K{i,k} * y{i,k+1} + g_ik - f{i,k};
                            for j = 1:i
                                r_ik = r_ik + Ca{i,j,k} * y{j,k+1};
                            end
                            r_ik = r_ik - Cg{i,k} * y{Ns, k};
                            
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
                        
                        d_ik = (-gamma(i) / hk) * y{end,k};
                        for j = 1:i
                            d_ik = d_ik + (alpha(i,j) / hk) * y{j,k+1};
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
                    display(sprintf('Iteration %d, Residual = %0.3g', nShooting, shootingRes));
                end

                if shootingRes > shootingTol
                    if isSymmetric
                        A = @ShootingNewton.MVPStoredLDL;
                        if verbose
                            dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Cg, L, D, P, S, Ns, Nt);
                            display(sprintf(' '));
                        else
                            [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, K, Ca, Cg, L, D, P, S, Ns, Nt);
                        end
                    else
                        A = @ShootingNewton.MVPStoredLU;
                        if verbose
                            dy = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Cg, L, U, P, Q, R, Ns, Nt);
                            display(sprintf(' '));
                        else
                            [dy,~,~,~] = gmres(A, r, maxGMRES, gmresTol, 1, [], [], r, Ca, Cg, L, U, P, Q, R, Ns, Nt);
                        end
                    end
                    y{end,1} = y{end,1} - dy;
                end
            end
        end
    end
    
    methods (Static)
        function mvp = MVPStoredLDL(z_s0, Ca, Cg, L, D, P, S, Ns, Nt)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                for i = 1:Ns         
                    rhs = Cg{i,k} * z_sk;
                    for j = 1:(i - 1);
                        rhs = rhs - Ca{i,j,k} * z{j};
                    end
                    z{i} = S{i,k} * (P{i,k} * (L{i,k}' \ (D{i,k} \ (L{i,k} \ (P{i,k}' * (S{i,k} * rhs))))));
                end
                z_sk = z{Ns};
            end
            mvp = z_sk - z_s0;
        end
        
        function mvp = MVPStoredLU(z_s0, Ca, Cg, L, U, P, Q, R, Ns, Nt)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                for i = 1:Ns         
                    r_ik = Cg{i,k} * z_sk;
                    for j = 1:(i - 1);
                        r_ik = r_ik - Ca{i,j,k} * z{j};
                    end
                    z{i} = Q{i,k} * (U{i,k} \ (L{i,k} \ (P{i,k} * (R{i,k} \ r_ik))));
                end
                z_sk = z{Ns};
            end
            mvp = z_sk - z_s0;
        end
        
        function mvp = MVP2(z_s0, Ca, Cg, J, Ns, Nt)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                for i = 1:Ns         
                    r_ik = Cg{i,k} * z_sk;
                    for j = 1:(i - 1);
                        r_ik = r_ik - Ca{i,j,k} * z{j};
                    end
                    z{i} = J{i,k} \ r_ik;
                end
                z_sk = z{Ns};
            end
            mvp = z_sk - z_s0;
        end
        
        function mvp = MVP0(z_s0, y, t, alpha, gamma, c, Ns, Nt, getMatrix)
            z_sk = z_s0;
            z    = cell(1, Ns);
            for k = 1:Nt
                hk = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    
                    C    = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                    r_ik = C * z_sk;
                    for j = 1:(i-1)
                        C    = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
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
        
        function [A,b,c] = getButcherTable(nStages)
            switch nStages
                case 1
                    A = 1;
                case 2
                    A = [1/3 0;
                         3/4 1/4];
                case 4
                  	A = [2  0 0 0;
                         3  3 0 0;
                         2 -2 4 0;
                         0  0 9 3]/12;
            end
            b = A(end, :);
            c = sum(A, 2);
        end
      
        function solverOut = configureSolver(varargin)
            if nargin > 0
                solverOut = ShootingNewton(varargin{:});
            end
        end
    end
end