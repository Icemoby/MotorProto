classdef (Sealed) Rect < Geometry2D
  	%Rect.m Creates an object representing a two dimensional rectangle
    %   G = Rect('PropertyName',propertyvalue,...) creates an object 
    %   representing a 2-Dimensional rectangle.
    %
    %   Example: Create a parameterized rectangle and change its properties.
    %       P = PARAMETER_LIST;
    %       P.new('len',1);
    %       P.new('wid',0.5);
    %       P.new('pos',rand(1,2));
    %       P.new('rot',pi/4);
    %       G = Rect('Length','len',...
    %                	'Width','wid',...
    %                   'Position','pos',...
    %                   'Rotation','rot');
    %       figure;subplot(1,2,1);axis equal;
    %       G.plot
    %
    %       P.edit('len',0.5);
    %       P.edit('wid',0.2);
    %       P.edit('pos',-rand(1,2));
    %       P.edit('rot',pi/3);
    %       G = rebuild(G);
    %       suplot(1,2,2);axis equal;
    %       G.plot;
    %
    % Rect methods:
    %   Rect defines no new methods
    %
    % Rect properties:
    %   Length - Length of the rectangle, oriented along the x-axis
    %   Width  - Width of the rectangle, oriented along the y-axis
    %
    % Rect inhertance:
    %   Rect inherts methods and properties from Geometry2D. See the help for
    %   Geometry2D for more information.
    %
    %   See also Geometry2D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    %% Public Properties
    properties (SetAccess = protected)
        Length
        Width
    end
    
    %% Public Methods
    methods
        %% Constructor Method
        function this = Rect(varargin)
            %% Set Defaults
            this             = this@Geometry2D;
            this.Domains     = {(1:4).'};
            
            %% Build Curves
            if nargin~=0               	
                nVarargin = numel(varargin);
                for i = 1:2:nVarargin
                    this.(varargin{i}) = varargin{i+1};
                end
                this = build(this);
            end
        end
        
        %% Parameterization Update Methods
        function this = refresh(this)
            warning('Rect:refresh','Replace "refresh" with "rebuild."');
            this = rebuild(this);
        end
        
        function this = rebuild(this)
            this          = this.rebuild@Geometry2D;
            this.Length   = rebuild(this.Length);
            this.Width    = rebuild(this.Width);
            this          = build(this);
        end
        
        function this = build(this)
            length   = this.Length;
            width    = this.Width;
            rotation = this.vRotation;
            position = this.vPosition;
            orient   = this.Orientation;
            
            switch this.Base
                case {'corner','Corner'}
                    points = [ 0      0;
                               length 0;
                               length width;
                               0      width];
                    
                    points = points*[ cos(rotation) sin(rotation);
                                     -sin(rotation) cos(rotation)];
                    
                    points(:,1) = points(:,1) + position(1);
                    points(:,2) = points(:,2) + position(2);
                case {'center','Center'}
                    points = [-length/2 -width/2;
                               length/2 -width/2;
                               length/2  width/2;
                              -length/2  width/2];
                    
                    points = points*[ cos(rotation) sin(rotation);
                                     -sin(rotation) cos(rotation)];
                    
                    points(:,1) = points(:,1) + position(1);
                    points(:,2) = points(:,2) + position(2);
                otherwise
                    error('Rect:build',...
                    'Invalid option %s for property Base in class rectangle',...
                    this.Base);
            end
            
            curveArray(1)	= Geometry1D.draw('Polyline',...
                                              'Points',points([1 2],:),...
                                              'Orientation',orient);
            curveArray(2) 	= Geometry1D.draw('Polyline',...
                                              'Points',points([2 3],:),...
                                              'Orientation',orient);
            curveArray(3)  	= Geometry1D.draw('Polyline',...
                                              'Points',points([3 4],:),...
                                              'Orientation',orient);
            curveArray(4)  	= Geometry1D.draw('Polyline',...
                                              'Points',points([4 1],:),...
                                              'Orientation',orient);
            this.Curves     = curveArray;
        end
    end
end