classdef (Sealed) Sector < Geometry2D
  	%Sector.m Creates an object representing a two dimensional annulus
    %   G = Sector('PropertyName',propertyvalue,...) creates an object 
    %   representing a 2-Dimensional annullus.
    %
    %   Example: Create a parameterized sector and change its properties.
    %       P = PARAMETER_LIST;
    %       P.new('rad',[0.5 1]);
    %       P.new('ang',pi/2);
    %       P.new('pos',rand(1,2));
    %       P.new('rot',pi/4);
    %       G = Sector('Radius','rad',...
    %                   'Angle','ang',...
    %                   'Position','pos',...
    %                   'Rotation','rot');
    %
    %       figure;subplot(1,2,1);axis equal;
    %       G.plot
    %
    %       P.edit('rad',[0.75 1]);
    %       P.edit('ang',pi/3);
    %       P.edit('pos',-rand(1,2));
    %       P.edit('rot',pi/3);
    %       G = rebuild(G);
    %       suplot(1,2,2);axis equal;
    %       G.plot;
    %
    % Sector methods:
    %   Sector defines no new methods
    %
    % Sector properties:
    %   Radius - An array containing the inner and outer radius of the annulus
    %   Angle  - The angle subtended by the annulus
    %
    % Sector inhertance:
    %   Sector inherts methods and properties from Geometry2D. See the help for
    %   Geometry2D for more information.
    %
    %   See also Geometry2D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    properties (SetAccess=protected)
        Radius
        Angle
    end
    
    methods
        function this = Sector(varargin)
            this.Domains     = {(1:4).'};
            if nargin ~= 0               
                for i = 1:2:nargin
                    this.(varargin{i}) = varargin{i+1};
                end
                this = build(this);
            end
        end
        
        function this = set.Radius(this,radiusIn)
            this.Radius = this.setProperty(radiusIn);
        end
        
        function this = set.Angle(this,angleIn)
            this.Angle = this.setProperty(angleIn);
        end
        
        function this = rebuild(this)
            this          = this.rebuild@Geometry2D;
            this.Radius   = rebuild(this.Radius);
            this.Angle    = rebuild(this.Angle);
            this          = build(this);
        end
        
    	function this = build(this)
            radius   = this.Radius.Value;
            angle    = this.Angle.Value;
            rotation = this.vRotation;
            position = this.vPosition;
            orient   = this.Orientation;
            base     = this.Base;
            switch base
                case 'center'
                    pointArray = [ 	radius(1)              0;
                                    radius(2)              0;
                                    radius(2)*cos(angle)  radius(2)*sin(angle);
                                    radius(1)*cos(angle)  radius(1)*sin(angle)];
                    
                    pointArray = pointArray*[	 cos(rotation) sin(rotation);
                                                -sin(rotation) cos(rotation)];
                    
                    pointArray(:,1) = pointArray(:,1)+position(1);
                    pointArray(:,2) = pointArray(:,2)+position(2);
                otherwise
                    error('Sector:build',...
                        'Invalid option for property Base in class Sector');
            end
            if abs(angle) < 2*pi*(1-eps)
                curveArray(1) = Geometry1D.draw('Polyline',...
                                	'Points',pointArray([1 2],:),...
                                    'Orientation',orient);
                curveArray(2) = Geometry1D.draw('Arc',...
                                	'Radius',radius(2),...
                                    'Angle',angle,...
                                    'Rotation',rotation,...
                                    'Position',position,...
                                    'Base',base,...
                                    'Orientation',orient);
                curveArray(3) = Geometry1D.draw('Polyline',...
                                    'Points',pointArray([3 4],:));
             	if radius(1) > 0 
                	curveArray(4) = Geometry1D.draw('Arc',...
                                    	'Radius',radius(1),...
                                    	'Angle',-angle,...
                                      	'Rotation',rotation+angle,...
                                       	'Position',position,...
                                       	'Base',base,...
                                      	'Orientation',orient);
                    this.Domains = {(1:4).'};
                else
                    this.Domains = {(1:3).'};
                end
            else
                curveArray(1) = Geometry1D.draw('Arc',...
                                                'Radius',radius(2),...
                                                'Angle',angle,...
                                                'Rotation',rotation,...
                                                'Position',position,...
                                                'Base',base,...
                                                'Orientation',orient);
                if radius(1) > 0
                    curveArray(2) = Geometry1D.draw('Arc', 'Radius'  , radius(1), 'Angle', -angle, 'Rotation'   , rotation+angle,...
                                                           'Position', position , 'Base' , base  , 'Orientation', orient);
                    this.Domains  = {(1:2).'};
                else
                    this.Domains = {1};
                end
                
            end
            this.Curves   = curveArray;
        end
    end
end