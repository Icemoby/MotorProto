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
    % See also MotorProto, 
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
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
        TimePoints = 0;
        Times      = 0;
        
        Verbose    = false;
    end
    
    properties (Dependent, SetAccess = private)
        Harmonics
    end
    
    properties (SetAccess = protected)
        X
        X_t
        
     	SimulationTime
        ConvergenceHistory
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