classdef Wire < handle & Parameterizable & matlab.mixin.Heterogeneous
    %Wire.m A concrete class representing insulated conductor models
    %   Wire objects
    %
    % Wire properties:
    %   ConductorMaterial - Material for the wire's conductor
    %   InsulatorMaterial - Material for the wire's insulation
    %
  	% Wire methods:
    %   build - Divides an input region into a number of conductors
    %
    % See also MotorProto, Model, Assembly, Stator, Slot
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%ConductorMaterial - Material for the wire's conductor
    %
    % See also Wire
    ConductorMaterial;
    
 	%InsulatorMaterial - Material for the wire's insulation
    %
    % See also Wire
    InsulatorMaterial;
%}
    properties
        ConductorMaterial = CopperExampleMaterial;
        InsulatorMaterial = Air;
    end
    
    methods
        %% Setters
        function this = set.ConductorMaterial(this, value)
            assert(isa(value, 'MaterialProperty'), 'MotorProto:Wire', 'ConductorMaterial must be a MaterialProperty object');
            this.ConductorMaterial = value;
        end
        
        function this = set.InsulatorMaterial(this, value)
            assert(isa(value, 'MaterialProperty'), 'MotorProto:Wire', 'InsulationMaterial must be a MaterialProperty object');
            this.InsulatorMaterial = value;
        end
    end
    
    methods (Abstract)
    	slotGeometry = build(slotShape, conductorBoundaries, nTurns, conductorDynamics, label)
    end
end