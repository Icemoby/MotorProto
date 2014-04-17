classdef SelfExcitedSynchronousRotor < PoleAndToothAssembly
    properties
        Slot
        OperatingMode = 'synchronous'
    end
    
    properties (Dependent)
        ConductorDynamics
        ConductorMaterial
        InsulatorMaterial

        WindingType
        Layers
        Turns
        
        FieldSlots
        TransformerSlots
        
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
     	function this = SelfExcitedSynchronousRotor(varargin)
            this = this@PoleAndToothAssembly(varargin{:});
            
            if isempty(this.Slot)
                this.Slot = Slot;
            end
            
        	newCircuit = Component.newComponent('FieldWoundTransformer', [this.Name,' Transformer']);
            this = addCircuit(this, newCircuit);
        end
        
        %% Getters
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
            value = this.Poles / 2;
        end
        
        function value = get.HasHalfWaveSymmetry(~)
        	value = true;
        end
        
        function value = get.GeometricSymmetries(this)
            value = this.Teeth;
        end
        
        function value = get.SpaceTimeSymmetries(this)
        	value = this.Poles / 2;
        end
        
        function value = get.AngularVelocity(this)
            switch this.OperatingMode
                case 'synchronous'
                    value = 4 * pi * this.ElectricalFrequency / this.Poles;
                case 'locked'
                    value = 0;
            end
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
        
        function value = get.FieldSlots(this)
            value = this.Circuits.FieldSlots;
        end
        
        function value = get.TransformerSlots(this)
            value = this.Circuits.TransformerSlots;
        end
        
        %% Setters
        function this = set.Slot(this, value)
            assert(isa(value, 'Slot'), 'MotorProto:Slot', 'SelfExcitedSynchronousRotor.Slot must be a Slot object');
            this.Slot = value;
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
        
        function this = set.FieldSlots(this, value)
            this.Circuits.FieldSlots = value;
        end
        
        function this = set.TransformerSlots(this, value)
            this.Circuits.TransformerSlots = value;
        end
        
        %% Others
        function [conductors, nonConductors, connectionMatrix] = buildPreProcessing(this)
            [conductors, nonConductors, connectionMatrix] = build(this.Slot, this.Name);
        end
        
        function this = buildPostProcessing(this, connectionMatrix, modeledFraction)
            nTurnsPerSlot           = this.Slot.Turns;
            nSlotsPerField          = this.FieldSlots / this.Poles;
            nSlotsPerTransformer    = this.TransformerSlots / this.Poles;
            nSlotsPerPole           = this.Teeth / this.Poles;
            nModeledPoles           = round(this.Poles * modeledFraction);
            
            connectionMatrices = cell(1, 2);
            connectionPolarity = cell(1, 2);
            
            for i = 1:nModeledPoles
                IField                = (1:nSlotsPerField) + (i-1) * nSlotsPerPole;
                subMatrix              = reshape(connectionMatrix(:,:,IField), [], nTurnsPerSlot * nSlotsPerField).';
                connectionMatrices{1} = [connectionMatrices{1}; subMatrix];
                connectionPolarity{1} = [connectionPolarity{1}; (-1)^(i+1) * ones(size(subMatrix))];
                
                ITransformer          = (1:nSlotsPerTransformer) + (i-1) * nSlotsPerPole + nSlotsPerField;
                subMatrix             = reshape(connectionMatrix(:,:,ITransformer), [], nTurnsPerSlot * nSlotsPerTransformer).';
                connectionMatrices{2} = [connectionMatrices{2}; subMatrix];
                connectionPolarity{2} = [connectionPolarity{2}; (-1)^(i+1) * ones(size(subMatrix))];
            end

            this.Circuits.ConnectionMatrices = connectionMatrices;
            this.Circuits.ConnectionPolarity = connectionPolarity;
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
            assemblyOut = SelfExcitedSynchronousRotor(varargin{:});
        end
    end
end