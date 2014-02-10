classdef PoleAndToothAssembly < RotatingMachineAssembly
    %PoleAndToothAssembly.m An abstract class representing assemblies defined by poles and teeth
    %   PoleAndToothAssembly objects are assemblies where the fundamental
    %   frequency of the spatial distribution and the number of fundamental
    %   geometric units do not coincide.
    %
    % PoleAndToothAssembly properties:
    %   Poles - The number of poles of the object
    %   Teeth - The number of teeth of the object
    %
    % PoleAndToothAssembly inherits properties and methods RotatingMachineAssembly.
    %
    % See the help files for RotatingMachineAssembly for more information.
    %
    % See also MotorProto, Model, RotatingMachineAssembly, Stator
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%Poles - The number of poles of the object
    %   The Poles property indicates the fundamental spatial frequency of the
    %   excitation and fields in the object. Specifically, half-wave symmetry is
    %   assumed and the fundamental frequency is Poles / 2.
    %
    % See also PoleAndToothAssembly
    Poles;
    
 	%Teeth - The number of teeth of the object
    %   The Teeth property indicates the number of fundamental geometric units
    %   required to represent the entire object. One geometric unit occupies an
    %   angle of 2*pi/Teeth.
    %
    % See also PoleAndToothAssembly
    Teeth;
%}
    properties
        Poles = PoleAndToothAssembly.setProperty(4);
        Teeth = PoleAndToothAssembly.setProperty(12);
    end

  	properties (Dependent)
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        GeometryFrequency
    end
    
    methods
        %% Constructor
     	function this = PoleAndToothAssembly(varargin)
            this = this@RotatingMachineAssembly(varargin{:});
        end
        
        %% Setters
        function set.Teeth(this,valueIn)
            this.Teeth = PoleAndToothAssembly.setProperty(valueIn);
        end       
        
        function set.Poles(this,valueIn)
            this.Poles = PoleAndToothAssembly.setProperty(valueIn);
        end

        %% Getters
        function value = get.GeometryFrequency(this)
            warning('Use GeometricSymmetries instead');
        	value = this.Teeth.Value;
        end

        function value = get.SolutionSpatialFrequency(this)
            warning('Use SpatialSymmetries instead');
        	value = this.Poles.Value  / 2;
        end
        
        function value = get.SolutionHalfWaveSymmetry(~)
            warning('Use HasHalfWaveSymmetry instead');
        	value = true;
        end
    end
    
   	methods (Static)
        function assemblyOut = newAssembly(varargin)
            assemblyOut = PoleAndToothAssembly(varargin{:});
        end
    end
end