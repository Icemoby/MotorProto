classdef SynchronousRotor < PoleAssembly
    %SynchronousRotor.m A concrete class representing Rotors intendid to operate at synchronous speed.
    %   SynchronousRotor objects represent assemblies intendid to operate at
    %   synchronous speed. A locked rotor simulation mode is also available by
    %   setting the OperatingMode property.
    %
    % SynchronousRotor properties:
    %   OperatingMode   - Sets one of two operating states of the object
    %   AngularVelocity - Speed of the rotor in radians per second.
    %
    % SynchronousRotor inherits properties and methods PoleAssembly.
    %
    % See the help files for PoleAssembly for more information.
    %
    % See also MotorProto, Model, PoleAssembly
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%OperatingMode - Sets one of two operating states of the object
    %   The OperatingMode property is used to toggle between one of two
    %   simulation configurations:
    %   
    %   {Synchronous} - The rotor is spinning in synchronicity with the 
    %                   fundamental harmonic. The angular
    %                   velocity is determined by the ElectricalFrequency and
    %                   Poles property.
    %
    %   Locked        - The rotor is locked and the angular velocity is zero.
    %
    % See also SynchronousRotor, AngularVelocity
    OperatingMode;
    
 	%AngularVelocity - Speed of the rotor in radians per second.
    %   The AngularVelocity is the speed of the rotor in radians per second. It
    %   varies depending on the OperatingMode property. The synchronous speed is
    %   defined as ElectricalFrequency * Poles / 2.
    %
    % See also SynchronousRotor, OperatingMode
    ConnectionType;
%}

    properties (Dependent)
        SolutionSpaceTimeSymmetry
        SolutionSpaceTimeCoefficients
        AngularVelocity
    end
    
    properties
        OperatingMode = 'synchronous';
    end
    
    methods
        %% Constructor
     	function this = SynchronousRotor(varargin)
            this = this@PoleAssembly(varargin{:});
        end
        
        %% Getters
        function value = get.SolutionSpaceTimeSymmetry(this)
            warning('MotorProto:Verbose', 'Use SpaceTimeSymmetries instead')
            value = [this.Poles, inf];
        end
        
        function value = get.SolutionSpaceTimeCoefficients(this)
            value = [-1, 1, 1];
        end
        
        function value = get.AngularVelocity(this)
            switch this.OperatingMode
                case 'synchronous'
                    value = 4 * pi * this.ElectricalFrequency / this.Poles;
                case 'locked'
                    value = 0;
            end
        end
        
        %% Setters
        function this = set.OperatingMode(this, value)
            switch lower(value)
                case {'synchronous', 'locked'}
                    this.OperatingMode = lower(value);
                otherwise
                    error('MotorProto:SynchronousRotor', 'Operating mode must be either  "synchronous" or "locked"');
            end
        end
    end
    
   	methods (Static)
        function assemblyOut = newAssembly(varargin)
            assemblyOut = SynchronousRotor(varargin{:});
        end
    end
end