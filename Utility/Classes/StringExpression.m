classdef StringExpression < Expression
    %StringExpression Stores an expression which depends on the PARAMETER_LIST
    %   S = StringExpression(D) creates a StringExpression object which stores a
    %   string that can be evaluated with the PARAMETER_LIST. The contents
    %   of D may be a string or a symbolic expression (requires the MATLAB
    %   Symbolic Math Toolbox).
    %
 	%   Example: StringExpression and ValueExpressions can be used
    %            interoperably, but ValueExpressions must be numeric.
    %       P = PARAMETER_LIST;
    %       P.new('myParam','rand(1)');
    %       S = StringExpression('myParam')
    %       V = ValueExpression(S.Value)
    %       P = P.rebuild
    %       S = S.rebuild
    %       V = V.rebuild
    %       V.Definition = S.Value;
    %       V
    %
    % StringExpression methods:
    %   rebuild     - Recalculates the value property
    %
    % StringExpression properties:
    %   Definition	- Holds the string which defines the expression
    %   Value       - Holds the most recent evaluation of the expression
    %
    %	See also ValueExpression, Parameterizable, PARAMETER_LIST, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
%{
properties
    %Definition - Holds the string which defines the expression
    %   Definition is a string which can be evaluated using the MATLAB eval
    %   function and the variables in the PARAMETER_LIST. It may be a scaler but
    %   can also contain function calls, vectors, and matrices.
    %
    % See also eval, Parameterizable, PARAMETER_LIST
    Definition;
    
    %Value - Holds the most recent evaluation of the expression
    %   Value contains the most recent evaluation of the Definition property
    %   with the PARAMETER_LIST. If Value is not up to date, the rebuild method
    %   may be called to re-evaluate the expression.
    %
    % See also rebuild, Parameterizable, PARAMETER_LIST.
    Value;
%}

    properties (SetAccess = protected)
        Definition
        Value
    end
    
    methods
        %% Constructor
        function this = StringExpression(definitionIn)
            if nargin~=0
                this.Definition = definitionIn;
            end
        end
        
        %% Set Methods
        function this = set.Definition(this,definitionIn)
            switch class(definitionIn)
                case {'sym'}
                    %% Convert Matrices
                    [nRows,mColumns]   = size(definitionIn);
                    charDefinition = '[';
                    for i = 1:nRows
                        charDefinition	= [charDefinition,...
                                          	char(definitionIn(i,1))];
                        for j = 2:mColumns
                            charDefinition  = [charDefinition,...
                                                ',',char(definitionIn(i,j))];
                        end
                        charDefinition	= [charDefinition,';'];
                    end
                    charDefinition 	= [charDefinition,']'];
                    this.Definition	= charDefinition;
                case {'char'}
                    this.Definition = definitionIn;
                otherwise
                    error('StringExpression:Definition',...
                            'Invalid input type %s for argument defIn',...
                            class(definitionIn));
            end
            this = rebuild(this);
        end
        
     	function this = rebuild(this)
            %rebuild - Recalculates the value property
            %   S1 = rebuild(S2) Creates a new StringExpression object S2 from 
            %   the StringExpression object S2. S1 and S2 have the same 
            %   definitions but S2's value property is set by re-evaluating the 
            %   Definition property using the PARAMETER_LIST's current state.
            %
            % See also StringExpression, PARAMETER_LIST, Parameterizable
            
            %% get the current paramList workspace
            PARAMETER_LIST.getWorkspace;
            
            %% get the number of String Expressions
            nExpressions = numel(this);
            
            for iExpression = 1:nExpressions;
                %% try evaluating the object definition
                try
                    newValue = eval(this(iExpression).Definition);
                catch ME
                    %% make sure all variables are defined in the PARAMETER_LIST
                    variableArray                = symvar(this.Definition);
                    nVariables                   = length(variableArray);
                    variablesExist(nVariables,1) = true;
                    
                    for iVariable = 1:nVariables
                        if ~exist( variableArray{iVariable}, 'var' )
                            variablesExist(iVariable) = false;
                        end
                    end
                    
                    %% check for errors
                    if any( ~variablesExist )
                        missingVariables    = variableArray(~variablesExist);
                        nMissing            = length(missingVariables);
                        errorMessage        = '"';
                        for iMissing = 1:nMissing
                            errorMessage = [errorMessage, missingVariables{iMissing}, '", '];
                        end
                        errorMessage(end-1:end) = [];
                        error('strExpression:refresh', 'Unknown parameters %s', errorMessage);
                    else
                        ME.rethrow
                    end
                end
                
                %% assign output
                this(iExpression).Value = newValue;
            end
        end
    end
end