classdef (HandleCompatible) Parameterizable
    %Parameterizable.m A mixin used in linking classes to the PARAMETER_LIST
    %   Users may define parameter's in the global PARAMETER_LIST in order to
    %   easily modify an existing design. The Parameterizable class provides a
    %   template for definining objects which depend on user defined parameters.
    %   Each parameterizable classs implements a method called rebuild which may
    %   be run to recalculate object properties after the user changes parameter
    %   values.
    %
    % Parameterizable methods:
    %   setProperty - Set a class properties as a variable expression
    %
    % See also StringExpression, ValueExpression, PARAMETER_LIST, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    %% Methods
    methods (Static,Access=protected)
        function expressionObjectOut = setProperty(definitionIn)
            %setProperty - Set a class properties as a variable expression
            %   P = setProperty('Expression') is intended for use within the 
            %   property set methods for Parameterizable subclasses. The 
            %   property P is then set as an object and which can update its 
            %   value from the user defined parameters in the PARAMETER_LIST.
            %
            % See also StringExpression, ValueExpression, PARAMETER_LIST
            
            switch class(definitionIn)
                case {'StringExpression','ValueExpression'}
                    expressionObjectOut = definitionIn;
                otherwise
                    if isnumeric(definitionIn)
                        expressionObjectOut = ValueExpression(definitionIn);
                    else
                        expressionObjectOut = StringExpression(definitionIn);
                    end
            end
        end
    end
    
    methods
        function this = rebuild(this)
            %% Rebuilds all properties which are defined as string expressions
            properties  = meta.class.fromName(class(this)).PropertyList;
            nProperties = numel(properties);
            for iProperty = 1:nProperties
                propertyName = properties(iProperty).Name;
                if isa(this.(propertyName),'StringExpression')
                    this.(propertyName) = rebuild(this.(propertyName));
                end
            end
        end
    end
end