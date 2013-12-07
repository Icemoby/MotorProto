classdef PoleAndToothAssembly < RotatingMachineAssembly  
    properties
        Poles = PoleAndToothAssembly.setProperty(2);
        Teeth = PoleAndToothAssembly.setProperty(6);
    end

  	properties (Dependent)
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        GeometryFrequency
    end
    
    methods
     	function this = PoleAndToothAssembly(varargin)
            this = this@RotatingMachineAssembly(varargin{:});
        end
        
        function set.Teeth(this,valueIn)
            this.Teeth = PoleAndToothAssembly.setProperty(valueIn);
        end       
        
        function set.Poles(this,valueIn)
            this.Poles = PoleAndToothAssembly.setProperty(valueIn);
        end

        function intOut = get.GeometryFrequency(this)
        	intOut = this.Teeth.Value;
        end

        function intOut = get.SolutionSpatialFrequency(this)
        	intOut = this.Poles.Value  / 2;
        end
        
        function boolOut = get.SolutionHalfWaveSymmetry(~)
        	boolOut = true;
        end
    end
    
   	methods (Static)
        function assemblyOut = newAssembly(varargin)
            assemblyOut = PoleAndToothAssembly(varargin{:});
        end
    end
end