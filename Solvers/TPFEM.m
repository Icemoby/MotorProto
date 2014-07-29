classdef TPFEM < Solver
    %TPFEM.m Time-Periodic FEM steady-state solver
    %
    % TPFEM properties:
    % 	NewtonTolerance       - Tolerance of iteration occuring at each time point
    % 	GMRESTolerance        - Tolerance of the initial condition correction
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
        NewtonTolerance     = 1e-6;
        GMRESTolerance      = 1e-6;
        MaxNewtonIterations = 100;
        MaxGMRESIterations  = 100;
        RungeKuttaStages    = 2;
        StorageLevel        = 0;
        SymmetricJacobian   = true;
    end
    
    methods
        %% Constructor
        function this = TPFEM(varargin)
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
            
            %% Get Algorithm Parameters
            Ns        = this.RungeKuttaStages;
            [A, ~, c] = this.getButcherTable(Ns);
            
            alpha = inv(A);
            gamma = sum(alpha,2);
            
            %% Initialize
            t          = model.getTimePoints(this.TimePoints);
            this.Times = t;
            Nt         = numel(t) - 1;

            r = cell(Ns, Nt);
            [r{:}] = deal(zeros(Nx, 1));
            if nargin == 3
                y = y0(:,1:Nt);
            else
                y = cell(Ns, Nt);
                [y{:}] = deal(zeros(Nx, 1));
            end
            
            %% Begin Simulation Timing
            if this.Verbose
                display(sprintf('TPFEM %d/%d\n',Ns,Nt));
            end
            
            tic;
            switch this.StorageLevel
                case 0
                    [y, r] = solve_level_0(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix);
                case 1
                    error('#TODO');
                case 2
                    [y, r] = solve_level_2(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix);
                case 3
                	[y, r] = solve_level_3(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix);
                otherwise
                    error('Storage level must be between 0 and 3');
            end
            this.SimulationTime = toc;
            
           	if this.Verbose
                display(sprintf('Simulation Time = %0.3g seconds\n', this.SimulationTime));
            end
            
            %% Post Processing
           	x   = cell(Nt,1);
            x_t = cell(Nt,1);
            for i = 1:Nt
                hk = t(i+1) - t(i);
                
                j = mod(i-2,Nt)+1;
                x{i} = y{end,j};
                
                k = mod(j-2,Nt)+1;
                x_t{i} = (- gamma(end) / hk) * y{end,k};
                for k = 1:Ns
                    x_t{i} = x_t{i} + (alpha(end,k) / hk) * y{k,j};
                end
            end
            x{Nt+1}   = x{1};
            x_t{Nt+1} = x_t{1};
            
            [x, x_t] = getMatrix.doPostProcessing(x, x_t);
            this.Y   = y;
            this.X   = x;
            this.X_t = x_t;
            solution = Solution(this);
        end
        
        function [y, r] = solve_level_0(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix)
            maxNewton = this.MaxNewtonIterations;
            maxGMRES  = this.MaxGMRESIterations;
            newtonTol = this.NewtonTolerance;
            gmresTol  = this.GMRESTolerance;
            verbose   = this.Verbose;
            
            nNewton = 0;
            newtonRes = 1;
            while newtonRes > newtonTol && nNewton < maxNewton
                nNewton = nNewton + 1;
                normf = 0;
                normr = 0;
                for k = 1:Nt
                    km1 = mod(k-2, Nt) + 1;
                    hk = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        K_ik = getMatrix.K(t_ik, h_ik);
                        r{i,k} = K_ik * y{i,k};
                        
                        f_ik = getMatrix.f(t_ik, h_ik);
                        r{i,k} = r{i,k} - f_ik;
                        
                        [~, g_ik] = getMatrix.G(t_ik, y{i,k}, h_ik);
                        r{i,k} = r{i,k} + g_ik;
                        
                        for j = 1:i
                            C      = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                            r{i,k} = r{i,k} + C * y{j,k};
                        end
                        C      = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                        r{i,k} = r{i,k} - C * y{Ns, km1};

                        normf = normf + norm(f_ik)^2;
                        normr = normr + norm(r{i,k})^2;
                    end
                end
                newtonRes = sqrt(normr / normf);

                if verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', nNewton, newtonRes));
                end
                
                if newtonRes > newtonTol
                    dy = reshape(r, [], 1);
                    dy = cell2mat(dy);
                    A = @TPFEM.MVP0;
                    
                    if verbose
                        dy = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, y, t, alpha, gamma, c, Ns, Nt, Nx, getMatrix);
                        display(sprintf(' '));
                    else
                        [dy,~,~,~] = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, y, t, alpha, gamma, c, Ns, Nt, Nx, getMatrix);
                    end
                    dy = TPFEM.PC0(dy, y, t, alpha, gamma, c, Ns, Nt, Nx, getMatrix);
                    dy = mat2cell(dy, Nx * ones(1, Nt*Ns));
                    dy = reshape(dy, Ns, Nt);
                    
                    for i = 1:Nt
                        for j = 1:Ns
                            y{j,i} = y{j,i} - dy{j,i};
                        end
                    end
                end
            end
        end
        
        function [y, r] = solve_level_2(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix)
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            isSymmetric = this.SymmetricJacobian;
            
            verbose = this.Verbose;
            
            K  = cell(Ns, Nt);
            f  = cell(Ns, Nt);
            Ca = cell(Ns, Ns, Nt);
            Cg = cell(Ns, Nt);
            J  = cell(Ns, Nt);
            
            for k = 1:Nt
                hk = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    
                    K{i,k}  = getMatrix.K(t_ik, h_ik);
                    Cg{i,k} = getMatrix.C(t_ik, hk / gamma(i), h_ik); 
                    f{i,k}  = getMatrix.f(t_ik, h_ik);
                    
                    for j = 1:i
                        Ca{i,j,k} = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                    end
                end
            end
                    
            nNewton = 0;
            newtonRes = 1;
            while newtonRes > newtonTol && nNewton < maxNewton
                nNewton = nNewton + 1;
                normf = 0;
                normr = 0;
                for k = 1:Nt
                    km1 = mod(k-2, Nt) + 1;
                    hk = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k}, h_ik);
                        
                        r{i,k} = K{i,k} * y{i,k} + g_ik - f{i,k};
                        for j = 1:i
                            r{i,k} = r{i,k} + Ca{i,j,k} * y{j,k};
                        end
                        r{i,k} = r{i,k} - Cg{i,k} * y{Ns, km1};
                        
                        J{i,k} = K{i,k} + Ca{i,i,k} + G_ik;
                        
                        normf = normf + norm(f{i,k})^2;
                        normr = normr + norm(r{i,k})^2;
                    end
                end
                newtonRes = sqrt(normr / normf);

                if verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', nNewton, newtonRes));
                end
                
                if newtonRes > newtonTol
                    dy = reshape(r,[],1);
                    dy = cell2mat(dy);
                    
                    A = @TPFEM.MVP2;
                    if verbose
                        dy = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, J, Ns, Nt, Nx);
                        display(sprintf(' '));
                    else
                        [dy,~,~,~] = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, J, Ns, Nt, Nx);
                    end
                    dy = TPFEM.PC2(dy, Ca, Cg, J, Ns, Nt, Nx);

                    dy = mat2cell(dy,Nx*ones(1,Nt*Ns));
                    dy = reshape(dy,Ns,Nt);
                    
                    for i = 1:Nt
                        for j = 1:Ns
                            y{j,i} = y{j,i} - dy{j,i};
                        end
                    end
                end
            end
        end
        
        function [y, r] = solve_level_3(this, y, r, t, alpha, gamma, c, Nt, Ns, Nx, getMatrix)
            maxNewton   = this.MaxNewtonIterations;
            maxGMRES    = this.MaxGMRESIterations;
            newtonTol   = this.NewtonTolerance;
            gmresTol    = this.GMRESTolerance;
            isSymmetric = this.SymmetricJacobian;
            
            verbose = this.Verbose;
            
            K  = cell(Ns, Nt);
            f  = cell(Ns, Nt);
            Ca = cell(Ns, Ns, Nt);
            Cg = cell(Ns, Nt);
            
            if isSymmetric 
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
                hk = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                    
                    K{i,k}  = getMatrix.K(t_ik, h_ik);
                    Cg{i,k} = getMatrix.C(t_ik, hk / gamma(i), h_ik); 
                    f{i,k}  = getMatrix.f(t_ik, h_ik);
                    
                    for j = 1:i
                        Ca{i,j,k} = getMatrix.C(t_ik, hk / alpha(i,j), h_ik);
                    end
                end
            end
                    
            nNewton = 0;
            newtonRes = 1;
            while newtonRes > newtonTol && nNewton < maxNewton
                nNewton = nNewton + 1;
                normf = 0;
                normr = 0;
                for k = 1:Nt
                    km1 = mod(k-2, Nt) + 1;
                    hk = t(k+1) - t(k);
                    for i = 1:Ns
                        t_ik = t(k) + c(i) * hk;
                        h_ik = hk / alpha(i,i);
                        
                        [G_ik, g_ik] = getMatrix.G(t_ik, y{i,k}, h_ik);
                        
                        r{i,k} = K{i,k} * y{i,k} + g_ik - f{i,k};
                        for j = 1:i
                            r{i,k} = r{i,k} + Ca{i,j,k} * y{j,k};
                        end
                        r{i,k} = r{i,k} - Cg{i,k} * y{Ns, km1};
                        
                        J_ik = K{i,k} + Ca{i,i,k} + G_ik;
                        if isSymmetric
                            [L{i,k}, D{i,k}, P{i,k}, S{i,k}] = ldl(J_ik);
                        else
                            [L{i,k}, U{i,k}, P{i,k}, Q{i,k}, R{i,k}] = lu(J_ik);
                        end
                        
                        normf = normf + norm(f{i,k})^2;
                        normr = normr + norm(r{i,k})^2;
                    end
                end
                newtonRes = sqrt(normr / normf);

                if verbose
                    display(sprintf('Iteration %d, Residual = %0.3g', nNewton, newtonRes));
                end
                
                if newtonRes > newtonTol
                    dy = reshape(r,[],1);
                    dy = cell2mat(dy);
                    
                    if isSymmetric
                        A = @TPFEM.MVPStoredLDL;
                        if verbose
                            dy = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, L, D, P, S, Ns, Nt, Nx);
                            display(sprintf(' '));
                        else
                            [dy,~,~,~] = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, L, D, P, S, Ns, Nt, Nx);
                        end
                        dy = TPFEM.PCStoredLDL(dy, Ca, Cg, L, D, P, S, Ns, Nt, Nx);
                    else
                        A = @TPFEM.MVPStoredLU;
                        if verbose
                            dy = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, L, U, P, Q, R, Ns, Nt, Nx);
                            display(sprintf(' '));
                        else
                            [dy,~,~,~] = gmres(A, dy, maxGMRES, gmresTol, 1, [], [], dy, Ca, Cg, L, U, P, Q, R, Ns, Nt, Nx);
                        end
                        dy = TPFEM.PCStoredLU(dy, Ca, Cg, L, D, P, S, Ns, Nt, Nx);
                    end
                    dy = mat2cell(dy,Nx*ones(1,Nt*Ns));
                    dy = reshape(dy,Ns,Nt);
                    
                    for i = 1:Nt
                        for j = 1:Ns
                            y{j,i} = y{j,i} - dy{j,i};
                        end
                    end
                end
            end
        end
    end
    
    methods (Static)
        %% Level 0
     	function y = MVP0(x, z, t, alpha, gamma, c, Ns, Nt, Nx, getMatrix)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            y = x;
            
            %% k = 1 Iteration
            k  = 1;
            hk = t(k+1) - t(k);
            for i = 1:Ns
                t_ik = t(k) + c(i) * hk;
                h_ik = hk / alpha(i,i);
                
                J_ik = getMatrix.K(t_ik, h_ik);
                J_ik = J_ik + getMatrix.C(t_ik, h_ik, h_ik);
                J_ik = J_ik + getMatrix.G(t_ik, z{i,k}, h_ik);
                
                x{i,k} = J_ik \ x{i,k};
                
                for j = (i+1):Ns
                    C      = getMatrix.C(t_ik, hk / alpha(j,i), hk / alpha(j,j));
                    x{j,k} = x{j,k} - C * x{i,k};
                end
            end
            
            %% k > 1 Iteration
            for k = 2:Nt
                hk = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                
                    C = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                    x{i,k} = x{i,k} + C * x{Ns, k-1};
                    
                    J_ik = getMatrix.K(t_ik, h_ik);
                    J_ik = J_ik + getMatrix.C(t_ik, h_ik, h_ik);
                    J_ik = J_ik + getMatrix.G(t_ik, z{i,k}, h_ik);
                    
                    x{i,k} = J_ik \ x{i,k};
                    
                    for j = (i+1):Ns
                        C      = getMatrix.C(t_ik, hk / alpha(j,i), hk / alpha(j,j));
                        x{j,k} = x{j,k} - C * x{i,k}; 
                  	end
                end
            end
            
            %% Matrix contribution from upper-right hand corner block
            k  = 1;
            hk = t(k+1) - t(k);
            for i = 1:Ns
                t_ik = t(k) + c(i) * hk;
                h_ik = hk / alpha(i,i);
                
              	C      = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                y{i,1} = y{i,1} - C * x{Ns, Nt};
            end
            
            y = reshape(y, Ns * Nt, 1);
            y = cell2mat(y);
        end
        
        function x = PC0(x, z, t, alpha, gamma, c, Ns, Nt, Nx, getMatrix)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            y = x;
            
            %% k = 1 Iteration
            k  = 1;
            hk = t(k+1) - t(k);
            for i = 1:Ns
                t_ik = t(k) + c(i) * hk;
                h_ik = hk / alpha(i,i);
                
                J_ik = getMatrix.K(t_ik, h_ik);
                J_ik = J_ik + getMatrix.C(t_ik, h_ik, h_ik);
                J_ik = J_ik + getMatrix.G(t_ik, z{i,k}, h_ik);
                
                x{i,k} = J_ik \ x{i,k};
                
                for j = (i+1):Ns
                    C      = getMatrix.C(t_ik, hk / alpha(j,i), hk / alpha(j,j));
                    x{j,k} = x{j,k} - C * x{i,k};
                end
            end
            
            %% k > 1 Iteration
            for k = 2:Nt
                hk = t(k+1) - t(k);
                for i = 1:Ns
                    t_ik = t(k) + c(i) * hk;
                    h_ik = hk / alpha(i,i);
                
                    C = getMatrix.C(t_ik, hk / gamma(i), h_ik);
                    x{i,k} = x{i,k} + C * x{Ns, k-1};
                    
                    J_ik = getMatrix.K(t_ik, h_ik);
                    J_ik = J_ik + getMatrix.C(t_ik, h_ik, h_ik);
                    J_ik = J_ik + getMatrix.G(t_ik, z{i,k}, h_ik);
                    
                    x{i,k} = J_ik \ x{i,k};
                    
                    for j = (i+1):Ns
                        C      = getMatrix.C(t_ik, hk / alpha(j,i), hk / alpha(j,j));
                        x{j,k} = x{j,k} - C * x{i,k}; 
                  	end
                end
            end
            
            x = reshape(x, Ns * Nt, 1);
            x = cell2mat(x);
        end
        
        %% Level 2
        function y = MVP2(x, Ca, Cg, J, Ns, Nt, Nx)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            y = x;
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:Ns
                x{j,k} = J{j,k} \ x{j,k};
                for i = (j+1):Ns
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:Nt
                for j = 1:Ns
                    x{j,k} = x{j,k} + Cg{j,k} * x{Ns,k-1};
                    x{j,k} = J{j,k} \ x{j,k};
                    for i = (j+1):Ns
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            %% Matrix contribution from upper-right hand corner block
            for j = 1:Ns
                y{j,1} = y{j,1} - Cg{j,1} * x{Ns,Nt};
            end
            
            y = reshape(y, Ns * Nt, 1);
            y = cell2mat(y);
        end
        
        function x = PC2(x, Ca, Cg, J, Ns, Nt, Nx)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:Ns
                x{j,k} = J{j,k} \ x{j,k};
                for i = (j+1):Ns
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:Nt
                for j = 1:Ns
                    x{j,k} = x{j,k} + Cg{j,k} * x{Ns,k-1};
                    x{j,k} = J{j,k} \ x{j,k};
                    for i = (j+1):Ns
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            x = reshape(x, Ns * Nt, 1);
            x = cell2mat(x);
        end
        
        %% Level 3
        function y = MVPStoredLDL(x, Ca, Cg, L, D, P, S, Ns, Nt, Nx)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            y = x;
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:Ns
                x{j,k} = S{j,k} * (P{j,k} * (L{j,k}' \ (D{j,k} \ (L{j,k} \ (P{j,k}' * (S{j,k} * x{j,k}))))));
                for i = (j+1):Ns
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:Nt
                for j = 1:Ns
                    x{j,k} = x{j,k} + Cg{j,k} * x{Ns,k-1};
                    x{j,k} = S{j,k} * (P{j,k} * (L{j,k}' \ (D{j,k} \ (L{j,k} \ (P{j,k}' * (S{j,k} * x{j,k}))))));
                    for i = (j+1):Ns
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            %% Matrix contribution from upper-right hand corner block
            for j = 1:Ns
                y{j,1} = y{j,1} - Cg{j,1} * x{Ns,Nt};
            end
            
            y = reshape(y, Ns * Nt, 1);
            y = cell2mat(y);
        end
        
        function x = PCStoredLDL(x, Ca, Cg, L, D, P, S, Ns, Nt, Nx)
            x = mat2cell(x, Nx * ones(1, Ns * Nt));
            x = reshape(x, Ns, Nt);
            
            %% k = 1 Iteration Setup
            k = 1;
            for j = 1:Ns
                x{j,k} = S{j,k} * (P{j,k} * (L{j,k}' \ (D{j,k} \ (L{j,k} \ (P{j,k}' * (S{j,k} * x{j,k}))))));
                for i = (j+1):Ns
                    x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k};
                end
            end
            
            %% k > 1 Iterations
            for k = 2:Nt
                for j = 1:Ns
                    x{j,k} = x{j,k} + Cg{j,k} * x{Ns,k-1};
                    x{j,k} = S{j,k} * (P{j,k} * (L{j,k}' \ (D{j,k} \ (L{j,k} \ (P{j,k}' * (S{j,k} * x{j,k}))))));
                    for i = (j+1):Ns
                        x{i,k} = x{i,k} - Ca{i,j,k} * x{j,k}; 
                  	end
                end
            end
            
            x = reshape(x, Ns * Nt, 1);
            x = cell2mat(x);
        end
        
        function y = MVPStoredLU(x, Ca, Cg, L, U, P, Q, R, Ns, Nt, Nx)
            % #TODO
        end
        
     	function y = PCStoredLU(x, Ca, Cg, L, U, P, Q, R, Ns, Nt, Nx)
            % #TODO
        end
        
        %%
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
                solverOut = TPFEM(varargin{:});
            end
        end
    end
end