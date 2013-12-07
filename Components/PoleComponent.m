classdef PoleComponent < RotatingMachineComponent    
    properties
        Poles = PoleComponent.setProperty(2)
    end

  	properties (Dependent)
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        GeometryFrequency
    end
    
    methods
     	function this = PoleComponent(varargin)
            this = this@RotatingMachineComponent(varargin{:});
        end

        function set.Poles(this,valueIn)
            this.Poles = PoleComponent.setProperty(valueIn);
        end
        
        function intOut = get.GeometryFrequency(this)
        	intOut = this.Poles.Value;
        end
        
        function intOut = get.SolutionSpatialFrequency(this)
        	intOut = this.Poles.Value  / 2;
        end
        
        function boolOut = get.SolutionHalfWaveSymmetry(~)
            boolOut = true;
        end
    end
    
   	methods (Static)
        function componentOut = newComponent(varargin)
            componentOut = PoleComponent(varargin{:});
        end
    end
end