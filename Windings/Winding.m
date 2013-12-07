classdef Winding < Parameterizable
    properties
        Phases = 3;
        Slots  = 6;
        Poles  = 2;
    end
    
    properties (Dependent, Abstract)
        WindingDiagram
    end
    
    methods
        %% Setters
        function this = set.Phases(this,value)
            this.Phases = this.setProperty(value);
        end
        
        function this = set.Slots(this,value)
            this.Slots = this.setProperty(value);
        end
        
        function this = set.Poles(this,value)
            this.Poles = this.setProperty(value);
        end
        
        %% Getters
        function value = get.Phases(this)
            value = this.Phases.Value;
        end
        
        function value = get.Slots(this)
            value = this.Slots.Value;
        end
        
        function value = get.Poles(this)
            value = this.Poles.Value;
        end
    end
end