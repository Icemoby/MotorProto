classdef Expression < matlab.mixin.Heterogeneous
    properties (Abstract,SetAccess = protected)
        Value
        Definition
    end
    
    methods (Abstract)
        this = rebuild(this)
    end
    
  	methods (Static,Access=protected,Sealed)
        function defaultObject = getDefaultScalarElement
            defaultObject      = ValueExpression(0);
        end
    end
end