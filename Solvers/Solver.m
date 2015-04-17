classdef Solver
    %Solver.m Base class for all solvers.
    %
    % Solver properties:
    %   Matrices           - MatrixFactory objects used by the solver
    %   TimePoints         - Target number of simulation time points
    %   Times              - The simulation times used by the solver
    %   Verbose            - Toggles verbose simulation progress reports
    %   X                  - The solution produced by the solver
    %   X_t                - Time derivative of the solver solution
    %   SimulationTime     - Total simulation time
    %   ConvergenceHistory - Relative residual history
    %
    % Solver methods:
    %   solve - Run the simulation
    %
    % See also MotorProto
    
%{
properties:
 	%TimePoints - Target number of simulation time points
 	%	The TimePoints property is an integer that sets the target number
 	%   of time points used in the simulation. This number is not adhered
 	%   to exactley but used as a guide to determine an appropriate number
 	%   of time points to use after considering harmonic content of the
 	%   solution and the effects of aliasing. In general, the number of
 	%   time points used will be greater than that specified in the
 	%   TimePoints, but not much larger.
    %
    % See also Solver, MatrixFactory/getTimePoints
    TimePoints;
%}

    properties (SetAccess = protected)
        Matrices
    end
    
    properties
        TimePoints
        Times
        
        StoreDecompositions = false;
        Verbose = false;
    end
    
    properties (Dependent, SetAccess = private)
        Harmonics
    end
    
    properties (SetAccess = protected)
        Y
        Y_t
        
        X
        X_t
        
     	SimulationTime
        DiscretizationError
    end
    
    methods (Static)        
        function solverOut = configureSolver(solverType,varargin)
            solverOut = eval(solverType);
            if isa(solverOut,'Solver')
                if nargin > 1;
                    solverOut = solverOut.configureSolver(varargin{:});
                end
            else
                error('MotorProto:Solver:invalidObjectType', '%s is not a recognized Solver subclass',solverType);
            end
        end
        
     	function [a,b,c,d,p,q,bu,be,pe] = getButcherTable(Ns)
            switch Ns
                case 1
                    a = [0 0;
                         0 1];
                    bu = [0;1];
                case 2
                    a = [0,             0,             0;
                         1 - 2^(1/2)/2, 1 - 2^(1/2)/2, 0;
                         2^(1/2)/4,     2^(1/2)/4,     1 - 2^(1/2)/2];
                    
                    bu = [sqrt(2)/2, -sqrt(2)/4;
                          sqrt(2)/2, -sqrt(2)/4;
                          1-sqrt(2), sqrt(2)/2];
                case 3            
                    c2 = roots([3 -18 18 -4]);
                    c2 = c2(3);
                    c2 = c2-((((c2-6)*c2+6)*c2)/(3*((c2-4)*c2+2))-4/(9*((c2-4)*c2+2)));
                    
                    a = zeros(4,4);
                    a(2,1) = c2 / 2;
                    a(2,2) = c2 / 2;
                    a(3,1) = (c2*(c2*(c2*(c2*(9*c2 - 30) + 38) - 20) + 4)) / (c2*(c2*(c2*(18*c2 - 48) + 56) - 32) + 8);
                    a(3,2) = (-c2*(c2*(c2*(3*c2 - 9) + 10) - 4)) / (c2*(c2*(c2*(9*c2 - 24) + 28) - 16) + 4);
                    a(3,3) = c2 / 2;
                    a(4,1) = (-c2*(c2*(c2*(9*c2 - 27) + 31) - 14) - 2) / (c2*(c2*(c2*(c2*(9*c2 - 54) + 102) - 84) + 24));
                    a(4,2) = (c2*(c2*(9*c2 - 30) + 34) - 12) / (c2*(c2*(12*c2 - 60) + 72) - 24);
                    a(4,3) = (-c2*(c2*(c2*(c2*(c2*(27*c2 - 108) + 198) - 208) + 132) - 48) - 8) / (c2*(c2*(c2*(c2*(c2*(36*c2 - 252) + 624) - 744) + 432) - 96));
                    a(4,4) = c2 / 2;
                    
                    bu = zeros(4,3);
                    bu(1,3) = -((3*c2^3 - 8*c2^2 + 10*c2 - 4))/(3*(3*c2^3 - 6*c2^2 + 4*c2)*(c2^2 - 4*c2 + 2));
                    bu(1,2) = ((9*c2^4 - 21*c2^3 + 15*c2^2 + 6*c2 - 6))/(3*(3*c2^3 - 6*c2^2 + 4*c2)*(c2^2 - 4*c2 + 2));
                    bu(1,1) = ((- 18*c2^4 + 51*c2^3 - 54*c2^2 + 18*c2))/(3*(3*c2^3 - 6*c2^2 + 4*c2)*(c2^2 - 4*c2 + 2));
                    
                    bu(2,3) = ((6*c2^2 - 20*c2 + 12))/(12*(c2 - 1)*(c2^2 - 4*c2 + 2));
                    bu(2,2) =  ((-9*c2^3 + 18*c2^2 + 6*c2 - 12))/(12*(c2 - 1)*(c2^2 - 4*c2 + 2));
                    bu(2,1) = ((18*c2^3 - 54*c2^2 + 48*c2 - 12))/(12*(c2 - 1)*(c2^2 - 4*c2 + 2));
                    
                    bu(3,3) = ((2*c2 - 4)*(3*c2^2 - 4*c2 + 2)^2)/(12*(c2^2 - 4*c2 + 2)*(- 3*c2^4 + 9*c2^3 - 10*c2^2 + 4*c2));
                    bu(3,2) = - ((3*c2^2 - 6)*(3*c2^2 - 4*c2 + 2)^2)/(12*(c2^2 - 4*c2 + 2)*(- 3*c2^4 + 9*c2^3 - 10*c2^2 + 4*c2));
                    bu(3,1) = - ((6*c2 - 6*c2^2)*(3*c2^2 - 4*c2 + 2)^2)/(12*(c2^2 - 4*c2 + 2)*(- 3*c2^4 + 9*c2^3 - 10*c2^2 + 4*c2));
                    
                    bu(4,3) = (2/(3*(c2^2 - 4*c2 + 2)));
                    bu(4,2) = (-(2*c2)/(c2^2 - 4*c2 + 2));
                    bu(4,1) = (c2^2/(c2^2 - 4*c2 + 2));
            end
            pe = size(bu,2);
            be = bu(:,pe).'*factorial(pe);
            
            b = a(end,:);
            c = sum(a,2);
            c = c(2:end);
            q = a(2:end,2:end) \ a(2:end,1);
            
            assert(abs(q(end)) < sqrt(eps));
            q(end) = 0;
            
            a = a(2:end,2:end);
            d = inv(a);
            p = sum(d,2);
        end
        
        function [z, z_t, s] = rkInterpolate(t, y, y_t, c, bu, s)
            [Rs, ~] = size(y);
            [m,n]   = size(bu);
            T       = t(end) - t(1);
            Nt      = length(t) - 1;
            
            Ns  = length(s) - 1;
            z   = cell(Rs, Ns);
            z_t = cell(Rs, Ns);
            for i = 1:Ns
                ds = s(i+1) - s(i);
                for j = 1:Rs
                    s_ij = s(i) + ds * c(j);
                    s_ij = mod(s_ij-T,T) + T;
                    s_ij = mod(s_ij-t(1),t(end)-t(1)) + t(1);
                    
                    k = find(t > s_ij, 1, 'first') - 1;
                    km1 = mod(k-2,Nt) + 1;
                    h_k = t(k+1)-t(k);
                    
                    u = (s_ij - t(k)) / h_k;

                    z{j,i} = y{Rs,km1};
                    z_t{j,i} = 0 * z{j,i};
                    for jj = 1:n
                        z{j,i} = z{j,i} + h_k * bu(1,jj) * y_t{Rs,km1} * u^jj;
                        z_t{j,i} = z_t{j,i} + bu(1,jj) * y_t{Rs,km1} * (jj*u^(jj-1)); 
                    end
                    
                    for ii = 2:m
                        for jj = 1:n
                            z{j,i} = z{j,i} + h_k * bu(ii,jj) * y_t{ii-1,k} * u^jj;
                            z_t{j,i} = z_t{j,i} + bu(ii,jj) * y_t{ii-1,k} * (jj*u^(jj-1)); 
                        end
                    end
                end
            end
        end
        
        function [ec, emax] = rkErrorCoefficients(t, y, y_t, be, pe, getMatrix)
            Nt = numel(t)-1;
            Ns = length(be) - 1;
            
            W = blkdiag(getMatrix.PostProcessing.SobolevA);
            W = W * blkdiag(getMatrix.PostProcessing.Reduced2Full);
            
            ec = zeros(1,Nt);
            Cn = 0;
            for k = 1:Nt
                km1 = mod(k-2,Nt) + 1;
                h_k = t(k+1) - t(k);
                
                ep = be(1) * W * y_t{end,km1};
                for i = 1:Ns
                   ep = ep + be(i+1) * W * y_t{i,k};
                end
                ep    = ep * h_k^(1-pe);
                ec(k) = norm(ep);
                
                Cn  = max(Cn, norm(W * y{end,k}));
            end
            ec = (ec ./ Cn) / factorial(pe);
            ec(1:(end/2)) = ec((end/2+1):end);
            ec = [ec(end),ec];
            h  = [t(end)-t(end-1), diff(t)];
            
            emax = max(ec.*h.^pe)^((pe+1)/pe);
            
%             figure;
%             plot(t, (ec.*(h.^(pe))).^((pe+1)/pe))
%             pause(1);
%             figure;
%             plot(t, ec)
%             pause(1);
        end
        
        function [s, tol] = rkRefine(t, tol, tol_max, ec, pe, rfact, init)
            tol_max = tol_max^(pe/(pe+1));
            
            if init
                h   = [t(end)-t(end-1),diff(t)];
                tol = mean(ec.*h.^pe);
                ngr = floor((-1/pe) * log(tol_max / tol) / log(rfact));
                tol = tol_max * (2^(pe*ngr));
            else
                tol = max(tol*2^(-pe), tol_max / 2);
            end
            
            Nt = length(t) - 1;
            T  = t(end) - t(1);
            
            if Nt < 12
                s = linspace(0,T,2*Nt+1);
            else
                %% Calculate new time-point function
                hc = (tol./ec).^(1/pe);
                I  = 2:(Nt/6+1);
                hc(I) = hc(I+3*Nt/6); %values in [T/2,T] are usually more well smoothed
                for i = 4:5
                    hc(I) = min(hc(I), hc(I+i*Nt/6));
                end

                for i = 1:5
                    hc(I+i*Nt/6) = hc(I);
                end
                hc(end) = hc(1);

                %% Calculate New Time Points
                t1 = t(1:(end-1));
                t2 = t(2:end);

                minax = [];
                s = [];
                imin = [];
                imax = [];
                for i = 2:(Nt-1)
                    j = i-1;
                    k = i+1;
                    if (hc(i) <= hc(j)) && (hc(i) < hc(k))
                        s     = cat(2,s,t(i));
                        minax = cat(2,minax,-1);
                        imin  = cat(2,imin,i);
                    elseif (hc(i) >= hc(j)) && (hc(i) > hc(k))
                        s     = cat(2,s,t(i));
                        minax = cat(2,minax,1);
                        imax  = cat(2,imax,i);
                    end
                end
                I = find(s > 0);
                if minax(I(1)) == -1
                    I(1) = [];
                end
                s = s(I);
                minax = minax(I);
                
                I = (s <= s(1) + T/6*(1+sqrt(eps)));
                s = s(I);
                minax = minax(I);
                assert(minax(end) == 1);
                
                m = length(s);
                
                for k = 1:2:m
                    if k < m %forward from max s(k) to min s(k+1)
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
                    end
                    
                    if k > 1 %backward from max s(k) to min s(k-1)
                        sb = [];
                        sk = s(k);
                        while sk > s(k-1)
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
                        sb = (sb - s(k)) * (s(k) - s(k-1)) / (s(k)-sb(end)) + s(k);
                        sb(end) = [];
                        s = cat(2,s,sb);
                    end
                end

                s = sort(s);
                s(end) = [];
                s = [s,s+T/6,s+T/3,s+T/2,s+2*T/3,s+5*T/6];
                s = mod(s,T);
                s = sort(s);
                s = [s,s(1)+T];
            end
        end
    end
    
    methods
        function value = get.Harmonics(this)
            if ~isempty(this.Matrices)
                value = this.Matrices.getHarmonics(this.Times);
            else
                hMax = floor(numel(this.Times) / 2);
                if hMax > 0
                    value = [0:hMax -hMax:1:-1];
                else
                    value = [];
                end
            end
        end
        
        function X = fft(this, x)
            t = this.Times;
            dt = diff(t);
            if all(dt-mean(dt)) < sqrt(eps) * mean(dt)
                X = fft(x,[],2) / (numel(t) - 1);
            else
                a = 2*pi*t / t(end);
                N = numel(a)-1;
                if mod(N,2) == 0
                    K = [0:(N/2) (1-N/2):-1];
                    D = exp(1i*K.'*a(1:end-1));
                    D(N/2+1,:) = cos(N/2*a(1:end-1));
                else
                    K = [0:((N-1)/2), ((1-N)/2):-1];
                    D = exp(1i*K.'*a(1:end-1));
                end
                X = x/D;
            end
    	end
    end
    
    methods (Sealed)
        function copyOut = copy(this)
            nThis   = numel(this);
            copyOut = this;
            for i = 1:nThis
                copyOut(i).Matrices = copy(this.Matrices);
            end
        end
    end
    
    methods (Abstract)
        solution = solve(this, model, x0);
    end
end