classdef (Sealed) Polygon2D < Geometry2D
  	%Polygon2D.m Creates an object representing a two dimensional polygon
    % G = Polygon2D('PropertyName',propertyvalue,...) creates an object
    % representing a two dimensional polygon.
    %
    %   Example: Create a parameterized polygon and change its properties.
    %       P = PARAMETER_LIST;
    %       P.new('pts',[0 0;1 0;0 1]);
    %       P.new('pos',[0 0]);
    %       P.new('rot',0);
    %
    %       G = Polygon2D('Points','pts','Rotation','rot','Position','pos')
    %       figure;subplot(1,2,1);axis equal;
    %       G.plot;
    %       
    %       P.edit('pts',[0 0;1 0;1 1;0 1]);
    %       P.edit('pos',rand(1,2));
    %       P.edit('rot',pi*rand);
    %       
    %       G = rebuild(G);
    %       G.PlotStyle = {'b'};
    %       subplot(1,2,2);axis equal;
    %       G.plot;
    %
    % Polygon2D methods:
    %   Polygon2D defines no new methods
    %
    % Polygon2D properties:
    %   Points - An n by 2 matrix representing the vertices of the polygon
    %
    % Polygon2D inhertance:
    %   Polygon2D inherts methods and properties from Geometry2D. See the help 
    %   for Geometry2D for more information.
    %
    % See also Geometry2D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    %% Public Properties
    properties (SetAccess=protected)
        Points
    end
    
    %% Public Methods
    methods
        %% Constructor Method
        function this = Polygon2D(varargin)
            this = this@Geometry2D;
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
            warning('Polygon2D:refresh','Change "refresh" to "update"');
            this.Points      = refresh(this.Points);
            this             = build(this);
        end
        
       	function this = rebuild(this)
            this        = this.rebuild@Geometry2D;
            this.Points = rebuild(this.Points);
            this        = build(this);
        end
        
        function this = build(this)
            %% get parameters
            points   = this.Points;
            position = this.vPosition;
            rotation = this.vRotation;
            orient   = this.Orientation;
            
            %% rotate
            rotationMatrix = [ cos(rotation) sin(rotation);
                              -sin(rotation) cos(rotation)];
            points         = points * rotationMatrix;

            %% shift
            points(:,1)    = points(:,1)+position(1);
            points(:,2)    = points(:,2)+position(2);
            
            %% make Curves
            [nPoints,~]     = size(points);
            curveArray   	= Polyline.empty(0,1);
            for i = 1:nPoints-1
                curveArray(i)	= Geometry1D.draw('Polyline', ...
                                                  'Points', points([i i+1],:),...
                                                  'Orientation',orient);
            end
            curveArray(nPoints)	= Geometry1D.draw('Polyline',   ...
                                                'Points',points([nPoints 1],:),...
                                                  'Orientation',orient);
            this.Curves  = curveArray;
            this         = sortCurves(this);
        end
    end
end