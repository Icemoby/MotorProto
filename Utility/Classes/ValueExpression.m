classdef ValueExpression < Expression
    %ValueExpression Stores a value which is independent of the PARAMETER_LIST
    %   V = ValueExpression(D) creates a ValueExpression object which stores a
    %   the value in D in such a way that the use of the ValueExpression class 
    %   is interchangable with the StringExpression class. The contents of D 
    %   must be numeric (ISNUMERIC(D) must return true).
    %
    %   Example 1: ValueExpressions and StringExpression can be used
    %              interoperably, but ValueExpressions must be numeric.
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
    %   Example 2: ValueExpressions must be numeric.
    %       P = PARAMETER_LIST;
    %       P.new('myParam','rand(2,2)');
    %       V = ValueExpression('myParam');	%this statement will not run
    %
    % ValueExpression properties
    %   Value - Holds the value of the expression
    %
    %	See also StringExpression, Parameterizable, PARAMETER_LIST, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
%{
properties:
    %Value - Holds the value of the expression
    %   This property holds the value of the expression. It must be numeric.
    %   This property may be used much in the same way that the StringExpression
    %   value property is used.
    %
    % See also Parameterizable, StringExpression
%}
    
    properties (SetAccess = protected)
        Value
        Definition
    end

    methods
        function this = ValueExpression(valueIn)
            if nargin~=0
            	this.Value      = valueIn;
                this.Definition = valueIn;
            end
        end
    end

    methods
        function this = rebuild(this)
            %NOOP
        end
    end
end