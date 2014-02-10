classdef AsynchronousRotor < PoleAndToothAssembly
    properties
        Slip
    end
    
   	properties (Dependent)        
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
     	function this = AsynchronousRotor(varargin)
            this = this@PoleAndToothAssembly(varargin{:});
        end
        
        %% Getters       	
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
            value = this.Poles.Value / 2;
        end
        
        function value = get.SolutionSpaceTimeSymmetry(this)
            warning('MotorProto:Verbose', 'Use SpaceTimeSymmetries instead')
            value = [this.Poles.Value, inf];
        end
        
        function value = get.SolutionSpaceTimeCoefficients(this)
            value = [-1, 1, 1];
        end
        
        function value = get.Slip(this)
            value = this.Slip.Value;
        end
        
        function value = get.AngularVelocity(this)
        	value = 2 * pi * this.ElectricalFrequency.Value / (this.Poles.Value / 2) * (1 - this.Slip);
        end
        
        %% Setters
        function this = set.Slip(this, value)
            this.Slip = AsynchronousRotor.setProperty(value);
        end
    end
    
   	methods (Static)
        function assemblyOut = newAssembly(varargin)
            assemblyOut = AsynchronousRotor(varargin{:});
        end
    end
end