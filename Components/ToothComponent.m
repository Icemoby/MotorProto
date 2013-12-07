classdef ToothComponent < RotatingMachineComponent    
    properties
        Teeth = ToothComponent.setProperty(2)
    end

  	properties (Dependent)
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        GeometryFrequency
    end
    
    methods
     	function this = ToothComponent(varargin)
            this = this@RotatingMachineComponent(varargin{:});
        end
        
        function set.Teeth(this,valueIn)
            this.Teeth = ToothComponent.setProperty(valueIn);
        end

        function intOut = get.GeometryFrequency(this)
        	intOut = this.Teeth.Value;
        end

        function intOut = get.SolutionSpatialFrequency(this)
        	intOut = this.GeometryFrequency;
        end
        
        function intOut = get.SolutionHalfWaveSymmetry(this)
        	intOut = false;
        end
    end
    
   	methods (Static)
        function componentOut = newComponent(varargin)
            componentOut = ToothComponent(varargin{:});
        end
    end
end