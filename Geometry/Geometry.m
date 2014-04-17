classdef (HandleCompatible) Geometry <  matlab.mixin.Heterogeneous
    %%Geometry Abstract interface class for geometry objects
    %   This class defines a common interface for all objects representing the
    %   geometry of a model. It cannot be instantiated.
    %
    % Geometry methods:
    %   Geometry implements no new methods.
    %
    % Geometry properties:
    %   Name        - The name of the object
    %   Description - A short description of the purpose of the object
    %
  	% Geometry inheritance:
    %   Geometry inherts methods and properties from Parameterizable. Refer to 
    %   the help for Parameterizable for more details.
    %
    % See also Geometry0D, Geometry1D, Geometry2D, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
%{
properties:
    %Name - The name of the object
    %   The name of the object is meant to be a short hand reference to the
    %   object. It should be unique.
    %
    % See also Geometry
    Name;
    
    %Description - A short description of the purpose of the object
    %   The description of the object should contain enough prose to adequently
    %   inform someone unfamiliar with the current simulation what the intent of
    %   the Geometry is.
    %
    % See also Geometry
    Description;
%}
    
    %% Public Properties
    properties
        Name
        Description
    end
    
 	properties (Abstract)
        PlotStyle
    end
    
    properties (Abstract,Constant)
        Dimension
    end
    
    properties 
        bbRadius
        bbCenter
        pClass
    end
    
    properties (SetAccess = protected)
        Axis
    end

  	properties (Dependent,SetAccess=protected)
        Position
        Rotation
    end
    
    properties (SetAccess = protected)
        vPosition       = [0 0];
        vRotation       = 0;
        fCachingEnabled = false;
    end
    
    methods
        function this = Geometry(varargin)
            if nargin ~=0
               	nVarargin = numel(varargin);
                for i = 1:2:nVarargin
                    this.(varargin{i}) = varargin{i+1};
                end
            end
            this.pClass = class(this);
            this.Axis   = LocalAxis;
        end
        
        function this = set.Name(this,nameIn)
            assert(ischar(nameIn),'Geometry:Name:set',...
                    'Name must be a string');
            this.Name = nameIn;
        end
        
        function this = set.Description(this,descriptionIn)
            assert(ischar(descriptionIn),'Geometry:Description:set',...
                    'Description must be a string');
            this.Description = descriptionIn;
        end
        
       	function this = set.Rotation(this,angleIn)
            this.Axis.Rotation = angleIn;
            this.vRotation     = this.Axis.Rotation;
            if this.fCachingEnabled
                this = updateCache(this);
            end
        end
        
        function angleOut = get.Rotation(this)
            angleOut = this.Axis.Rotation;
        end
        
        function this = set.Position(this,xyIn)
            this.Axis.Position = xyIn;
            this.vPosition     = this.Axis.Position;
            if this.fCachingEnabled
                this = updateCache(this);
            end
        end
        
        function positionOut = get.Position(this)
            positionOut = this.Axis.Position;
        end
    end
    
    methods (Abstract)
        plot(obj)
        wireframe(obj)
    end
    
    methods (Abstract,Access=protected)
        this = updateCache(this)
    end
    
    methods
        function this = enableCaching(this)
            [this.fCachingEnabled] = deal(true);
        end
        
        function this = disableCaching(this)
            [this.fCachingEnabled] = deal(false);
        end
    end
    
    methods (Static, Abstract)
        geometryOut = draw(typeIn,varargin)
    end
    
  	methods (Static,Access=protected,Sealed)
        function defaultObject = getDefaultScalarElement
            defaultObject = Polyline;
        end
    end
    
    methods (Sealed)        
        function [I,J]   = boundingBallsOverlap(this1,this2)
            [X1,Y1,R1] = getBoundingBall(this1);
            if nargin == 1
                distanceToCenter = sqrt( bsxfun(@minus,X1,X1.').^2 ...
                                        +bsxfun(@minus,Y1,Y1.').^2);
                sumOfRadii       = bsxfun(@plus,R1,R1.');
                lessThanRadius   = bsxfun(@lt,distanceToCenter,...
                                              sumOfRadii*(1+sqrt(eps)));
                [I,J]            = find(tril(lessThanRadius,-1));
            else
                [X2,Y2,R2]       = getBoundingBall(this2);
                distanceToCenter = sqrt( bsxfun(@minus,X1,X2.').^2 ...
                                        +bsxfun(@minus,Y1,Y2.').^2);
                sumOfRadii       = bsxfun(@plus,R1,R2.');
                lessThanRadius   = bsxfun(@lt,distanceToCenter,...
                                              sumOfRadii*(1+sqrt(eps)));
                [I,J]            = find(lessThanRadius);
            end
        end
        
        function In      = inBoundingBall(this,xIn,yIn)
            [X,Y,R]          = getBoundingBall(this);
            distanceToCenter = sqrt( bsxfun(@minus,X,xIn.').^2 ...
                                    +bsxfun(@minus,Y,yIn.').^2);
            lessThanRadius   = bsxfun(@lt,distanceToCenter,R*(1+sqrt(eps)));
            In               = any(lessThanRadius,2);
        end
        
        function [X,Y,R] = getBoundingBall(this)
            n   = numel(this);
            Ix  = 1:2:(2*n);
            Iy  = 2:2:(2*n);             
            bbC = [this.bbCenter];
            X   = bbC(Ix).';
            Y   = bbC(Iy).';
            R   = [this.bbRadius].';
        end
        
        function classes = getElementClasses(this)
         	classes = {this.pClass};
        end
    end
end