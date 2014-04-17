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
        Slot
        ConnectionType = ConnectionTypes.Wye;
    end
    
    properties (Dependent)
        Phases
        SourceType        
        ConductorDynamics
        ConductorMaterial
        InsulatorMaterial

        WindingType
        Layers
        Turns
        
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
            
            if isempty(this.Slot)
                this.Slot = Slot;
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
            switch this.WindingType
                case WindingTypes.Distributed
                    value = this.Poles / 2;
                case WindingTypes.Concentrated
                    value = gcd(this.Poles,this.Teeth);
                    if (this.Layers == 1) && (value > 1)
                        value = value / 2;
                    end
                otherwise
                    error('No Implementation');
            end
        end
        
        function value = get.HasHalfWaveSymmetry(this)
            switch this.WindingType
                case WindingTypes.Distributed
                    value = true;
                case WindingTypes.Concentrated
                    g         = gcd(this.Poles,this.Teeth);
                    evenTeeth = (mod(this.Teeth / g, 2) == 0);
                    oddPoles  = (mod(this.Poles / g, 2) == 1);
                    if evenTeeth && oddPoles
                        value = true;
                    else
                        value = false;
                    end
                otherwise
                    error('No Implementation');
            end
        end
        
        function value = get.GeometricSymmetries(this)
            value = this.Teeth;
        end
        
        function value = get.SpaceTimeSymmetries(this)
            switch this.WindingType
                case WindingTypes.Distributed
                    value = this.Poles / 2 * this.Phases;
                case WindingTypes.Concentrated
                    error('No Implementation');
                otherwise
                    error('No Implementation');
            end
        end
        
        function value = get.SolutionSpaceTimeSymmetry(this)
            warning('MotorProto:Verbose','Use SpaceTimeSymmetries instead');
            value = [this.Teeth / 3, 3];
        end
        
        function value = get.SolutionSpaceTimeCoefficients(this)
            warning('MotorProto:Verbose', 'Combine this property and solutionspacetimesymmetry into a single property or multiple non-array properties');
            value = [1, 6, 2];
        end
        
        function value = get.AngularVelocity(~)
            value = 0;
        end
       
        function value = get.WindingType(this)
            value = this.Slot.WindingType;
        end
       
        function value = get.Layers(this)
            value = this.Slot.Layers;
        end
       
        function value = get.Turns(this)
            value = this.Slot.Turns;
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
               
        function this = set.WindingType(this, value)
            this.Slot.WindingType = value;
        end
               
        function this = set.Layers(this, value)
            this.Slot.Layers = value;
        end
               
        function this = set.Turns(this, value)
            this.Slot.Turns = value;
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
          	nTeeth = this.Teeth;
            nPoles = this.Poles;
            switch this.WindingType
                case WindingTypes.Distributed
                    assert(mod(nTeeth / nPoles, 3) == 0, 'MotorProto:StatorComponent', 'The number of teeth per pole must be an integer multiple of 3. The current value is %d', nTeeth/nPoles);
                case WindingTypes.Concentrated
                    assert(mod(nTeeth , 3) == 0,'MotorProto:StatorComponent', 'The number of teeth must be an integer multiple of 3. The current value is %d', nTeeth);
                otherwise
                    error('MotorProto:StatorComponent', 'Unknown WindingType %s', char(stator.WindingType));
            end
          	%% Perform Actions
            [conductors, nonConductors, connectionMatrix] = build(this.Slot, this.Name);
        end
        
        function this = buildPostProcessing(this, connectionMatrix, modeledFraction)
            nPhases = this.Phases;
            
            connectionMatrices = cell(1, nPhases);
            connectionPolarity = cell(1, nPhases);
            switch this.WindingType
                case WindingTypes.Distributed
                    nTurnsPerSlot  = this.Slot.Turns;
                    nSlotsPerPhase = this.Teeth / this.Poles / nPhases;
            
                    for i = 1:nPhases
                        I                     = (1:nSlotsPerPhase) + (i-1) * nSlotsPerPhase;
                        connectionMatrices{i} = reshape(connectionMatrix(:,:,I), [], nTurnsPerSlot * nSlotsPerPhase).';
                        if mod(i,2) == 1
                            connectionPolarity{i} =   ones(size(connectionMatrices{i}));
                        else
                            connectionPolarity{i} = - ones(size(connectionMatrices{i}));
                        end
                    end
                case WindingTypes.Concentrated
                    W = generateConcentratedWindingLayout(this.Poles, this.Teeth, this.Layers);
                    
                    if this.HasHalfWaveSymmetry && (mod(this.Poles * modeledFraction,2) == 1)
                        I = W(1:(end/2));
                        J = W((end/2+1):end);
                        J = (J == mod(I+2,6)+1);
                        
                        assert(all(J),'Winding layout does not result in half-wave symmetry');
                        W = I;
                        
                    end
                    
                    W = [W;mod(W+nPhases-1,2*nPhases)+1];
                    W = reshape(W,1,[]).';
                    
                    for i = 1:nPhases
                        if mod(i,2) == 1
                            Ip = (W == i);
                            In = (W == (i+nPhases));
                        else
                            Ip = (W == (i+nPhases));
                            In = (W == i);
                        end
                        
                        connectionMatrices{i} = cat(3,connectionMatrix(:,:,Ip),connectionMatrix(:,:,In));
                        connectionMatrices{i} = reshape(connectionMatrices{i},[],numel(connectionMatrices{i})).';
                        
                        connectionPolarity{i} = cat(3,1+0*connectionMatrix(:,:,Ip),-1+0*connectionMatrix(:,:,In));
                        connectionPolarity{i} = reshape(connectionPolarity{i},[],numel(connectionPolarity{i})).';
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