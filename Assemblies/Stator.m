classdef Stator < PoleAndToothAssembly
    %Stator.m A concrete class representing Stators with distributed windings
    %   Stator objects are assemblies representing integer stators with distrubuted
    %   windings having an integer number of slots per pole per phase.
    %
    % Stator properties:
    %   Slot           - An object representing the conductor layout and shape
    %                    of a single slot of the stator.
    %   ConnectionType - Indicates how the windings of the Stator are connected
    %   SourceType     - Indicates how the windings of the Stator are excited
    %
    % Stator inherits properties and methods PoleAndToothAssembly.
    %
    % Stator provides and interface for properties of Slot.
    %
    % See the help files for PoleAndToothAssembly and Slot for more information.
    %
    % See also MotorProto, Model, PoleAndToothAssembly, Slot, ConnectionTypes
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%Slot - An object representing the conductor layout and shape of a single slot of the stator.
    %   A Stator object's Slot property is a Slot object. The slot object
    %   specifies the shape of the slot and controls how the conductors are
    %   layed out within that area.
    %
    %   See the help for Slot for more information.
    %
    % See also Stator, Slot
    Slot;
    
 	%ConnectionType - Indicates how the windings of the Stator are connected
    %   The ConnectionType property indicates the configuration of the Stator
    %   windings, e.g. Wye/Delta. The property can be set using a string or
    %   through the ConnectionTypes enumeration object.
    %
    %   See the help for ConnectionTypes for more information.
    %
    % See also Stator, ConnectionTypes
    ConnectionType;
    
 	%SourceType - Indicates how the windings of the Stator are excited
    %   The SourceType property indicates the type of source connected to the
    %   stator windings, e.g. Current/Voltage. The property can be set using a
    %   string or through the SourceTypes enumeration object.
    %
    %   See the help for SourceTypes for more information.
    %
    % See also Stator, SourceTypes
    SourceTypes;
%}  
    properties
        Slot = Slot
        ConnectionType = ConnectionTypes.Wye;
    end
    
    properties (Dependent)
        Phases
        SourceType        
        ConductorDynamics
        ConductorMaterial
        InsulatorMaterial

        %% New Symmetry Properties
        SpatialSymmetries
        HasHalfWaveSymmetry
        SpaceTimeSymmetries
        GeometricSymmetries
        
        %% Old Symmetry Properties
        SolutionSpaceTimeSymmetry
        SolutionSpaceTimeCoefficients
        AngularVelocity
    end
    
    methods
        %% Constructor
     	function this = Stator(varargin)
            this = this@PoleAndToothAssembly(varargin{:});
            
            if isempty(this.Sources)
                this.SourceType = SourceTypes.CurrentSource;
            end
        end
        
        %% Getters
        function value = get.Phases(this)
            value = this.Sources.Phases;
        end
        
        function value = get.SourceType(this)
            value = this.Sources.Type;
        end
        
        function value = get.ConductorDynamics(this)
            value = this.Slot.ConductorDynamics;
        end
        
        function value = get.ConductorMaterial(this)
            value = this.Slot.ConductorMaterial;
        end
        
    	function value = get.InsulatorMaterial(this)
            value = this.Slot.InsulatorMaterial;
        end
        
       	function value = get.SpatialSymmetries(this)
            value = this.Poles.Value / 2;
        end
        
        function value = get.HasHalfWaveSymmetry(~)
            value = true;
        end
        
        function value = get.GeometricSymmetries(this)
            value = this.Teeth.Value;
        end
        
        function value = get.SpaceTimeSymmetries(this)
            value = this.Poles.Value / 2 * this.Phases.Value;
        end
        
        function value = get.SolutionSpaceTimeSymmetry(this)
            warning('MotorProto:Verbose','Use SpaceTimeSymmetries instead');
            value = [this.Teeth.Value / 3, 3];
        end
        
        function value = get.SolutionSpaceTimeCoefficients(this)
            warning('MotorProto:Verbose', 'Combine this property and solutionspacetimesymmetry into a single property or multiple non-array properties');
%             value = [1, 6, 6];
            value = [1, 6, 2];
        end
        
        function value = get.AngularVelocity(~)
            value = 0;
        end
        
        %% Setters
        function this = set.Phases(this, value)
            this.Sources.Phases = value;
        end
        
        function this = set.SourceType(this, value)
            if isa(value, 'SourceTypes')
                newType = value;
            elseif ischar(value)
                newType = SourceTypes.(value);
            else
                newType = SourceTypes(value);
            end
            
            oldSource = this.Sources;
            if isempty(oldSource)
                newSource   = Component.newComponent(char(newType), [this.Name,' Source']);
                this        = addSource(this, newSource);
            elseif ~(oldSource.Type == newType)
                newSource   = Component.newComponent(char(newType), this.Sources);
                this        = removeComponent(this, oldSource.Name);
                this        = addSource(this, newSource);
            end
        end
        
        function this = set.Slot(this, value)
            assert(isa(value, 'Slot'), 'MotorProto:Slot', 'Stator.Slot must be a Slot object');
            this.Slot = value;
        end
               
        function this = set.ConnectionType(this, value)
            if isa(value, 'ConnectionTypes')
                this.ConnectionType = type;
            elseif ischar(value)
                this.ConnectionType = ConnectionTypes.(value);
            else
                this.ConnectionType = ConnectionTypes(value);
            end
        end
        
        function this = set.ConductorDynamics(this, value)
            this.Slot.ConductorDynamics = value;
        end
        
        function this = set.ConductorMaterial(this, value)
            this.Slot.ConductorMaterial = value;
        end
        
        function this = set.InsulatorMaterial(this, value)
            this.Slot.InsulatorMaterial = value;
        end
        
        %% Others
        function [conductors, nonConductors, connectionMatrix] = buildPreProcessing(this)
            %% Check Configuration
          	nTeeth = this.Teeth.Value;
            nPoles = this.Poles.Value;
            assert(mod(nTeeth / nPoles, 3) == 0, 'MotorProto:StatorComponent', ...
                                               	 'The number of teeth per pole must be an integer multiple of 3. The current value is %d', nPoles/nTeeth);
            %% Perform Actions
            [conductors, nonConductors, connectionMatrix] = build(this.Slot, this.Name);
        end
        
        function this = buildPostProcessing(this, connectionMatrix)
            nTurnsPerSlot      = this.Slot.Turns;
            nPhases            = this.Phases.Value;
            nSlotsPerPhase     = this.Teeth.Value / this.Poles.Value / nPhases;
            
            connectionMatrices = cell(1, nPhases);
            connectionPolarity = cell(1, nPhases);
            for i = 1:nPhases
                I                     = (1:nSlotsPerPhase) + (i-1) * nSlotsPerPhase;
                connectionMatrices{i} = reshape(connectionMatrix(:,:,I), [], nTurnsPerSlot * nSlotsPerPhase).';
                if mod(i,2) == 1
                    connectionPolarity{i} =   ones(size(connectionMatrices{i}));
                else
                    connectionPolarity{i} = - ones(size(connectionMatrices{i}));
                end
            end
            this.Sources.ConnectionMatrices = connectionMatrices;
            this.Sources.ConnectionPolarity = connectionPolarity;
        end
        
        function previewElement(this)
            domainHull = makeElementDomainHull(this);
            
            if ~isempty(this.InputRegions)
                regionGeometry = this.InputRegions.Geometry;
                domainHull     = domainHull - regionGeometry;
            else
                regionGeometry = [];
            end
            
            if ~isempty(this.Slot.Shape)
                domainHull   = domainHull - this.Slot.Shape;
                turns        = this.buildPreProcessing;
                turnGeometry = [turns.Geometry];
            else
                turnGeometry = [];
            end
            
            plot([domainHull, turnGeometry, regionGeometry]);
        end
    end
    
    methods (Static)
        function assemblyOut = newAssembly(varargin)
            assemblyOut = Stator(varargin{:});
        end
    end
end