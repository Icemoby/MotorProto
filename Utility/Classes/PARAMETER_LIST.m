classdef (Sealed) PARAMETER_LIST < Parameter
    %% PARAMETER_LIST Global parameter list for parameterized simulations
    %	PARAMETER_LIST is used to define user input parameters on which
    %	various simulation objects may depend.
    %
    %   Example 1: Define and edit a parameter
    %       P = PARAMETER_LIST;
    %       P.new('myParam',1)
    %       P.edit('myParam',2)
    %
    %   Example 2: Define two parameters, where the second depends on the first
   	%       P = PARAMETER_LIST;
    %       P.new('myAngle',pi/2)
    %       P.new('myParam','myAngle*2')
    %       P.edit('myAngle',pi/3)
    %       P.Parameters(2)
    %
    %   Example 3: Define a parameter with a function call, and re-evaluate the
    %              function using the rebuild method
    %       P = PARAMETER_LIST;
    %       P.new('randParam','rand(1,2)')
    %       P.Parameters
    %       P.rebuild;
    %       P.Parameters
    %
   	%   Example 4: Define two different variables which refrence the unique
   	%              PARAMETER_LIST object.
   	%       P1 = PARAMETER_LIST;
    %       P1.new('myParam1',[-1,0;0,-1])
    %       P2 = PARAMETER_LIST;
    %       P1 == P2
    %       P2.new('myParam2','[0,1;1,0]*myParam1');
    %       P1.Parameters(1) == P2.Parameters(1)
    %       P1.Parameters(2) == P2.Parameters(2)
    %
    % PARAMETER_LIST methods:
    %   sort         	- Sorts parameters and checks for circular definitions
    %	new            	- Creates a new parameter
    %   edit          	- Changes the definition of an existing parameter
    %   cut             - Removes an existing parameter from the list
    %   rebuild      	- Updates parameter values from their definitions
    %   getWorkspace	- Assigns parameters to workspace with current values
    %   getSymbolic  	- Assigns parameters to workspace as symbolic variables
    %   
    % PARAMETER_LIST properties:
    %   Parameters      - An array of user defined parameters
    %
    % See also Parameter, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
%{
properties
	%Parameters - An array of user defined parameters.
 	%   Parameters is an array of Parameter objects, corresponding to the
 	%   user defined inputs. The order of appearance of the parameters in
  	%   should not be relied upon and may be different from the order in
 	%   which they have been defined. The sort method will reorder the
  	%   list so that parameters which depend on other members of the list
 	%   appear later in the list.
    %
    % See also PARAMETER_LIST
    Parameters;
end
%}
    
    %% Private Properties
    properties (SetAccess = private)
        Parameters = Parameter.empty(1,0);
    end
    
    %% Events
    events
        PARAMETER_LIST_UPDATE
    end
    
    %% Public Methods
    methods
        %% Singleton Constructor Method
        function THIS = PARAMETER_LIST
            persistent UNIQUE_INSTANCE
            if isempty(UNIQUE_INSTANCE)
                THIS.Name = 'Global Parameter Definitions';
                UNIQUE_INSTANCE = THIS;
            else
                THIS = UNIQUE_INSTANCE;
            end
        end

        %% Sort parameters and check for recursive parameter definitions
        function THIS = sort(THIS)
            %sort - Sorts parameters array and checks for circular definitions
            %   sort(P) Changes the order of apperance of the parameters in
            %   the parameter array P so that each parameter depends only on the
            %   parameters which precede it in the array. A warning will be
            %   issued if circular parameter definitions are detected.
            %
            % See also PARAMETER_LIST, Parameter
            
            %% get parameter names
            parameterArray  = THIS.Parameters;
            nameArray       = {parameterArray.Name};
            nParameters  	= numel(nameArray);
            
            %% find depdencies
            isDependentOn       = zeros(nParameters,nParameters);
            for iParameter      = 1:nParameters
                variableArray   = parameterArray(iParameter).Variables;
                nVariables      = numel(variableArray);
                for iVariable   = 1:nVariables
                    dependentVariables = strcmp(nameArray,...
                                                variableArray{iVariable});
                    isDependentOn(dependentVariables,iParameter) = 1;
                end
            end
            
            %% check for circular definitions
            %   The lu decomposition results in a permuation array which orders 
            %   the parameters so that the independent parameters are first, 
            %   parameters which are dependent on those are next, and so on. 
            %   The decomposition will also be upper triangular if there are no 
            %   recursive definitions, and is used for error checking.
            
            [~,~,permutationArray]  = lu(isDependentOn,'vector');
            isDependentOn           = isDependentOn(permutationArray,...
                                                    permutationArray);
            
            recursiveDefinitionsFound = any(any(tril(isDependentOn)));
            if recursiveDefinitionsFound
                warning('PARAMETER_LIST:sort',...
                        'Recursive parameter definitions found');
            end
            
            %% assign output
            THIS.Parameters         = parameterArray(permutationArray);
        end
        
        %% Disable concatentation (should never need to concatentate singleton)
        function vertcat(varargin)
            error( 'PARAM_LIST:vertcat',...
                'Concatentation of the global PARAMETER_LIST is not allowed.');
        end
        
        function horzcat(varargin)
            error(  'PARAM_LIST:vertcat',...
                'Concatentation of the global PARAMETER_LIST is not allowed.');
        end
    end
    
    %% Static Methods
    methods (Static)
        %% Add a new parameter to PARAMETER_LIST
        function parameterOut = new(nameIn,definitionIn)
            %new - Adds a new user defined parameter to the parameter list.
            %   P.new(N,D) adds a new user defined parameter named N with 
            %   definition D to the parameter list P. The parameter name should 
            %   be a character array. The parameter expression may be a 
            %   numerical value or character string containing an expression 
            %   for which the MATLAB eval function can operate on. If the 
            %   parameter expression contains the names of other parameters, 
            %   those parameters will be used when the expression is evaluated.
            %
            % See also eval, PARAMETER_LIST, Parameter, Parameterizable
            
            %% Check to see if the parameter exists
            THIS       	= PARAMETER_LIST;
            nameArray	= {THIS.Parameters.Name};
            iName       = find(strcmp(nameArray,nameIn));
            
            assert(isempty(iName),...
                    'PARAMETER_LIST:new',...
                    'A parameter named "%s" already exists',...
                    nameArray{iName});
            
            %% Construct Parameter
            THIS                  	= PARAMETER_LIST;
            parameterOut            = Parameter(nameIn,definitionIn);
            THIS.Parameters(end+1)	= parameterOut;
            
            %% Create callback for updates
            addlistener(parameterOut, 'Definition', 'PostSet', @THIS.rebuild);
            
            %% Assign symbolic variables to the base workspace and caller
            THIS.rebuild;
            assignin('base',    nameIn, sym(nameIn));
            assignin('caller',  nameIn, sym(nameIn));
        end
        
        %% Edit an existing parameter
        function parameterOut = edit(nameIn,definitionIn)
            %edit - Changes the definition of an existing parameter
            %   P.edit(N,D) will replace the existing definition of the 
            %   parameter named N with the definition D. A warning will be 
            %   issued if the new definition causes a circular definition.
            %
            % See also PARAMETER_LIST, Parameter, Parameterizable
            
            %% Check to see if the parameter exists
            THIS       	= PARAMETER_LIST;
            nameArray	= {THIS.Parameters.Name};
            iName       = find(strcmp(nameArray,nameIn));
            
            
            assert(~isempty(iName),...
                    'PARAMETER_LIST:edit',...
                    'No parameter named %s found',...
                    nameArray{iName});
            assert(numel(iName)==1,...
                    'PARAMETER_LIST:edit',...
                    'The parameter named %s is multiply defined',...
                    nameArray{iName(end)});

            %% Assign the new expression
            parameterOut                = THIS.Parameters(iName);
            parameterOut.Definition     = definitionIn;
            
            %% Check for recursive parameter definitions
          	THIS = sort(THIS);
            
            %% Refresh parameters and assign symbolic variables
            THIS.rebuild;
            assignin('base',    parameterOut.Name, sym(parameterOut.Name));
            assignin('caller',  parameterOut.Name, sym(parameterOut.Name));
        end
        
        %% Delete an existing parameter
        function cut(nameIn)
            %%cut - Removes an existing parameter from the list
            %   P.cut(N) will remove a parameter with name N from the parameter
            %   list if it exists. A warning will be issued if it does not.
            %   This method does nothing to ensure that no other objects depend
            %   on the parameter which is being removed. If a parameter is
            %   removed and objects still dependent on it, an error will be
            %   issued the next time the removed parameter is needed.
            %
            % See also PARAMETER_LIST, Parameter, Parameterizable
            
            %% Find parameter
            THIS            = PARAMETER_LIST;
            nameArray       = {THIS.Parameters.Name};
            iName           = find(strcmp(nameArray,nameIn));
            parameterExists = ~isempty(iName);
            
            %% Delete if it exists
            if parameterExists
                THIS.Parameters(iName) = [];
            else
                warning('PARAMETER_LIST:cut',...
                     	'No parameter named %s found.',...
                      	nameIn);
            end
            
            %% Recalculate the parameters
            THIS.rebuild;
        end
        
        %% Callback for execution when a parameter is changed, or manual refresh
        function THIS = refresh(~,~)
            warning('PARAMETER_LIST:refresh',...
                        'Replace "refresh" with "rebuild."');
            THIS = rebuild(PARAMETER_LIST);
        end
        
        function THIS = rebuild(~,~)
         	%rebuild - Updates parameter values from their definitions
            %   P.rebuild recalculates the value property of each parameter 
            %   object in the parameter list P from their definitions. This 
            %   method also sorts the parameters in order of dependency.
            %
            % See also sort PARAMETER_LIST, Parameter, Parameterizable
            
            %% setup
            THIS            = PARAMETER_LIST;
            THIS            = sort(THIS);
            parameterArray  = THIS.Parameters;
            nParameters     = numel(parameterArray);
            defArray        = cell(nParameters,1);
            valueArray      = cell(nParameters,1);
            varArray        = cell(nParameters,1);
            
            %% refresh each parameter
            for iParameter  = 1:nParameters
                parameterArray(iParameter).rebuild;
                defArray{iParameter} = parameterArray(iParameter).Definition;
                valueArray{iParameter} = parameterArray(iParameter).Value;
                varArray{iParameter} = parameterArray(iParameter).Variables;
            end
            
            %% assign cell arrays to the PARAMETER_LIST
            THIS.Definition = defArray;
            THIS.Value      = valueArray;
            THIS.Variables  = varArray;
            
            %% send notifcation that update is complete
            notify(THIS,'PARAMETER_LIST_UPDATE');
        end
        
        %% Assign parameter variables values to calling workspace
        function getWorkspace
            %getWorkspace - Assigns evaluated parameters to workspace
            %   P.getWorkspace creates variables in the calling functions 
            %   workspace from the parameter list P with names corresponding to 
            %   the names of the  parameters and values corresponding to the 
            %   current values of the parameters. If it is uncertain whether or
            %   not the parameter values are up to date, use the rebuild method
            %   before calling getWorkspace.
            %
            % See also rebuild, PARAMETER_LIST, Parameter, Parameterizable
            
            THIS            = PARAMETER_LIST;
            parameterArray  = THIS.Parameters;
            nParameters     = numel(parameterArray);
            for iParameter  = 1:nParameters
                assignin('caller',...
                            parameterArray(iParameter).Name,...
                            parameterArray(iParameter).Value);
            end
        end
        
        %% Assign parameters as symbolic variables to calling workspace
        function getSymbolic
           	%getSymbolic - Assigns symbolic parameters to workspace
            %   P.getSymbolic creates symbolic variables in the calling 
            %   function's workspace from the parameter list P with names 
            %   corresponding to the names of the parameters and symbolic 
            %   expressions which correspond to the parameter definitions. This
            %   method requires the MATLAB Symbolic Toolbox.
            %
            % See also PARAMETER_LIST, Parameter, Parameterizable
            
            THIS            = PARAMETER_LIST;
            parameterArray  = THIS.Parameters;
            nParameters     = numel(parameterArray);
            for iParameter  = 1:nParameters
                assignin('caller',....
                            parameterArray(iParameter).Name,...
                            sym(parameterArray(iParameter).Name));
            end
        end
    end
end