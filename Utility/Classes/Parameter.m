classdef Parameter < matlab.mixin.Copyable
    %Parameter.m A single user definable parameter
    %   Parameter objects cannot be created by themselves but are instantiated
    %   through the global PARAMETER_LIST. The Name property should be a string.
    %   The definition property may be set as either a numeric, string, or
    %   symbolic variable (requires the MATLAB Symbolic Toolbox).
    %
    % Parameter methods:
    %   rebuild     - Recalcultes the parameter's value property
    %
    % Parameter properties:
    %   Name        - String defining how the parameter should be referenced
    %   Definition  - The defining expression of the parameter
    %   Value       - The particular value of the evaluated definition
    %   Variables   - A list of variables on which the parameter depends
    %
    % See also PARAMETER_LIST, MOTOR_PROTO
    
 	%   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
%{
properties
    %Name - String defining how the parameter should be referenced
    %   Name defines the unique identifier for a particular parameter. This is
    %   the value that should be used to refer to the parameter in other
    %   parameter definitions and in parameterizable object properties
    %
    % See also Parameterizable
    Name;
    
    %Definition - The defining expression of the parameter
    %   Definition defines how the Value property of the Parameter should be
    %   calculated. Definition will typically be a string which can be evaluated
    %   using the MATLAB builtin function eval. In general, Definition can be
    %   assigned using any variable which is convertible to a character array.
    %
    % See also eval, Parameter, PARAMETER_LIST, MotorProto
    Definition;
    
    %Value - The particular value of the evaluated definition
    %   Value contains a variable which is created by evaluating the Parameter's
    %   Definition property against the current PARAMETER_LIST.
    %
    % See also Parameter, PARAMETER_LIST, MotorProto
    Value;
    
    %Variables - A list of variables on which the parameter depends
    %   Variables is a cell array containing strings representing the variables
    %   in the Definition property as recognized the MATLAB builtin symvar. For
    %   the Value property of the Parameter to be evaluated, these strings must
    %   correspond to Parameters which exist in the PARAMETER_LIST.
    %
    % See also symvar, Parameter, PARAMETER_LIST, MotorProto
    Variables;
%}
    
    %% Protected Properties
    properties (SetObservable,SetAccess = protected)
        Name
        Definition
    end
    
    %% Private Properties
    properties (SetAccess = protected)
        Value
        Variables
    end
    
    %% Protected Methods
    methods (Access = protected)
        %% Constructor method
        function this = Parameter(nameIn,definitionIn)
            if nargin ~= 0
                if nargin == 1
                    definitionIn = '';
                end
                this.Name       = nameIn;
                this.Definition = definitionIn;
            end
        end
    end
    
    %% Public Methods
    methods
        function rebuild(this)
            %%rebuild - Recalcultes the parameter's value property
            %   This method uses the current values of all the parameters
            %   defined in the parameter to list to update the value of this
            %   parameter. Each entry in the Variable property must correspond 
            %   to a Parameter in the PARAMETER_LIST in order for this method to
            %   successfuly run. 
            %
            % See also Parameter, PARAMETER_LIST, MotorProto
            
            if isnumeric(this.Definition)
                this.Value      = this.Definition;
                this.Variables  = '';
            else
                this.Variables  = symvar(this.Definition);
                PARAMETER_LIST.getWorkspace;
                this.Value      = eval(this.Definition);
            end
        end
        
        function set.Definition(this,definitionIn)
            switch class(this)
                case{'Parameter'}
                    if isnumeric(definitionIn)
                        this.Definition = definitionIn;
                    else
                        this.Definition = char(definitionIn);
                    end
                    this.rebuild;
                otherwise
                    this.Definition = definitionIn;
            end
        end
    end
end