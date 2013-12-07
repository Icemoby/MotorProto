classdef LocalAxis < Parameterizable
    properties (SetAccess = protected)
        InitialPosition = LocalAxis.setProperty([0 0]);
        InitialAngle    = LocalAxis.setProperty(0);
        RotationAngle   = LocalAxis.setProperty(0);
        RotationAxis    = LocalAxis.setProperty([0 0]);
    end
    
    properties (Dependent)
        Rotation
        Position
    end  
    
    properties
        cRotation = 0
        cPosition = [0 0]
    end
    
    properties (Dependent)
        Angle
    end
    
    methods
        %% Constructor Method
        function this = LocalAxis(varargin)
            if nargin~=0
                nVarargin = numel(varargin);
                for i = 1:2:nVarargin
                    this.(varargin{i}) = varargin{i+1};
                end
                this.Position  = this.calculatePosition;
                this.Rotation  = this.calculateRotation;
            end
        end
        
        %% Set Methods
        function this = set.InitialAngle(this,angleIn)
            this.InitialAngle = this.setProperty(angleIn);
        end
        
        function this = set.InitialPosition(this,positionIn)
            this.InitialPosition = this.setProperty(positionIn);
        end
        
        function this = set.Rotation(this,rotationAngle)
        	this.InitialAngle  = LocalAxis.setProperty(rotationAngle);
            this.RotationAngle = LocalAxis.setProperty(0);
         	this.cRotation     = calculateRotation(this);
        end
        
        function angleOut = get.Rotation(this)
            angleOut = this.cRotation;
        end
        
        function this = set.Position(this,rotationAxis)
        	this.InitialPosition = LocalAxis.setProperty(rotationAxis);
            this.RotationAxis    = LocalAxis.setProperty([NaN NaN]);
            this.cPosition       = calculatePosition(this);
        end
        
        function coordOut = get.Position(this)
            coordOut = this.cPosition;
        end
        
        function angleOut = get.Angle(this)
            angleOut = this.Rotation;
        end
        
        %% Cache Methods
        function angleOut = calculateRotation(this)
            angleOut = this.InitialAngle.Value + sum([this.RotationAngle.Value]);
            angleOut = mod(angleOut+pi,2*pi)-pi;
        end
        
        function positionOut = calculatePosition(this)
            positionOut   = this.InitialPosition.Value;
            rotationAngle = [this.RotationAngle.Value];
            nRotations    = numel(rotationAngle);
            rotationAxis  = reshape([this.RotationAxis.Value],[],nRotations).';
            for iRotation = 1:nRotations
                currentAxis    = rotationAxis(iRotation,:);
                if any(isnan(currentAxis))
                    currentAxis = positionOut;
                end
                currentAngle   = rotationAngle(iRotation);
                rotationMatrix = [  cos(currentAngle) sin(currentAngle);
                                   -sin(currentAngle) cos(currentAngle)];
                positionOut    =   currentAxis...
                                 +(positionOut-currentAxis)*rotationMatrix;
            end
        end
        
        %% Rotation Method
        function varargout = rotate(varargin)
            this    = [varargin{1:(end-2)}];
            angleIn = varargin{end-1};
            axisIn  = varargin{end};
            
            nThis    = numel(this);
            thisSize = size(this);
            
            nAngles             = numel(angleIn);
            newAngle(1,nAngles) = LocalAxis.setProperty(angleIn(end));
            if nAngles == 1
                newAngle = repmat(newAngle,thisSize);
            else
                for iAngle = 1:(nAngles - 1)
                    newAngle(1,iAngle) = LocalAxis.setProperty(angleIn(iAngle));
                end
            end
            angleIn = [newAngle.Value];
            newAngle = num2cell(newAngle);
            
            nAxis            = numel(axisIn) / 2;
            newAxis(1,nAxis) = LocalAxis.setProperty(axisIn(end,:));
            if nAxis == 1
                newAxis = repmat(newAxis,thisSize);
            else
                for iAngle = 1:(nAxis - 1)
                    newAxis(1,iAngle) = LocalAxis.setProperty(axisIn(iAngle,:));
                end
            end
            axisIn  = [newAxis.Value];
            axisIn  = reshape(axisIn,2,[]).';
            newAxis = num2cell(newAxis);
            
            %%
            newRotationAngle     = {this.RotationAngle};
            newRotationAngle     = cellfun(@(x,y)(vertcat(x,y)),...
                                            newRotationAngle,...
                                            newAngle,...
                                            'UniformOutput',false);
            [this.RotationAngle] = deal(newRotationAngle{:});            
                
            newRotationAxis     = {this.RotationAxis};
            newRotationAxis     = cellfun(@(x,y)(vertcat(x,y)),...
                                            newRotationAxis,...
                                            newAxis,...
                                            'UniformOutput',false);
            [this.RotationAxis] = deal(newRotationAxis{:});
            
         	%% Calculate the new angle
            newRotation      = [this.cRotation];
            newRotation      = newRotation + angleIn;
            newRotation      = mod(newRotation+pi,2*pi)-pi;
            newRotation      = num2cell(newRotation);
            [this.cRotation] = deal(newRotation{:});
            
            %% Calculate the new position
            newPosition            = [this.cPosition];
            newPosition            = reshape(newPosition,2,[]).';
            
            keepPosition           = any(isnan(axisIn),2);
            axisIn(keepPosition,:) = newPosition(keepPosition,:);
            newPosition            = newPosition - axisIn;
            newPosition            = [sum(newPosition .* [cos(angleIn.'),-sin(angleIn.')],2),...
                                      sum(newPosition .* [sin(angleIn.'), cos(angleIn.')],2)];
            newPosition            = newPosition + axisIn;
            newPosition            = mat2cell(newPosition,ones(nThis,1),2).';
            [this.cPosition]       = deal(newPosition{:});
            
            if nargout == 1
                varargout = {this};
            elseif nargout == nThis
                varargout = num2cell(this); 
            else
                varargout = [num2cell(this), newRotation, newPosition];
            end
        end
        
        %% Parameterization Update Methods
        function this = rebuild(this)
            this.InitialPosition = rebuild(this.InitialPosition);
            this.InitialAngle    = rebuild(this.InitialAngle);
            this.RotationAngle   = rebuild(this.RotationAngle);
            this.RotationAxis    = rebuild(this.RotationAxis);
            this.Position        = this.calculatePosition;
            this.Rotation        = this.calculateRotation;
        end
    end
end