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
        
     	function [a,b,c,d,p,q,bu,bh,bth] = getButcherTable(Ns)
            switch Ns
                case 1
                    a = [0 0;
                         0 1];
                         
                    bh  = [1,0];
                    bth = [1,0];
                case 2
                    a = [0,0,0;
                         1 - 2^(1/2)/2, 1 - 2^(1/2)/2, 0;
                         2^(1/2)/4,2^(1/2)/4, 1 - 2^(1/2)/2];
                    
                    bu  = [sqrt(2)/2,-sqrt(2)/4;
                           sqrt(2)/2,-sqrt(2)/4;
                           1-sqrt(2),sqrt(2)/2];
                           
                    bh  = [1,0,0];
                    bth = [1,0,0];
                case 3            
                    c2 = roots([3 -18 18 -4]);
                    c2 = c2(2);
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
                    
                    bh  = [ -((c2 - 1)^2*(3*c2^2 - 4*c2 + 2))/(c2*(3*c2^2 - 6*c2 + 4)*(c2^2 - 4*c2 + 2)), ((3*c2^2)/4 - 1/2)/((c2 - 1)*(c2^2 - 4*c2 + 2)) + 3/4, ((3*c2^2 - 4*c2 + 2)^2*(c2^2 - 2*c2 + 2))/(4*(c2^2 - 4*c2 + 2)*(- 3*c2^4 + 9*c2^3 - 10*c2^2 + 4*c2)), (c2*(c2 - 2))/(c2^2 - 4*c2 + 2)];
                    bth = [ ((3*c2 - 2)*(c2^2 - 2*c2 + 2))/(c2*(3*c2^2 - 6*c2 + 4)*(c2^2 - 4*c2 + 2)), - 1/(2*(c2 - 1)) - (c2 - 2)/(c2^2 - 4*c2 + 2), ((c2 - 2)*(3*c2^2 - 4*c2 + 2)^2)/(2*c2*(c2 - 1)*(3*c2^2 - 6*c2 + 4)*(c2^2 - 4*c2 + 2)), 1 - 2/(c2^2 - 4*c2 + 2)];
            end
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
        
        function ei = rkError(t, y, y_t, bh, W)
            Nt = length(t) - 1;
            Ns = length(bh) - 1;
           	ei = zeros(1, Nt);
            for k = 1:Nt
                km1 = mod(k-2, Nt) + 1;
                h_k = t(k+1) - t(k);
                yh  = y{end,km1} + h_k*bh(1)*y_t{end,km1};
                for i = 1:Ns
                    yh  = yh  + h_k*bh(i+1)*y_t{i,k};
                end
                
                ei(k) = sqrt((W*yh).'*(W*yh));
                yh    = yh - y{end,k};
                ei(k) = sqrt((W*yh).'*(W*yh)) / ei(k);
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