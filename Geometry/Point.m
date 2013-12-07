classdef Point < Geometry0D
    %%Point An object representing a point or set of points
    %   Point objects represent sets of points in the problem domain. They may
    %   be used to fix certain node locations in the resulting finite element
    %   mesh.
    %
    %   Example: Create a parameterized point
    %       P = PARAMETER_LIST;
    %       P.new('xParam','rand(1)');
    %       P.new('yParam','rand(1)');
    %
    %       G = Point('X','xParam','Y','yParam','PlotStyle',{'x','r'});
    %       figure;subplot(1,2,1);
    %       G.plot;
    %
    %       P.edit('xParam','-rand(1)');
    %       P.edit('yParam','-rand(1)');
    %       G = rebuild(G);
    %       subplot(1,2,2);
    %       G.plot;
    %
    % Point methods:
    %   rebuild - Recalculates the object's properties
    %
    % Point properties:
    %   X - An array of X coordinates
    %   Y - An array of Y coordinates
    %
    % Point inheritance:
    %   Point inherits all methods and properties of Geometry0D. Refer to the
    %   help for Geometry0D for more details.
    %
    % See also Geometry0D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $

%{
properties:
    %X - The X coordinate of the point
    %   This property represent the X coordinate (or set of coordinates) for the
    %	point object. It is parameterizable.
    %
    % See also Parameterizable
    X;
    
    %Y - The Y coordinate of the point
    %   This property represent the Y coordinate (or set of coordinates) for the
    %	point object. It is parameterizable.
    %
    % See also Parameterizable
    Y;
%}        
    %% Public Properties
    properties
        X
        Y
    end
    
    %% Public Methods
    methods
        %% Constructor
        function this = Point(varargin)
            %call superclass constructor
            this = this@Geometry0D(varargin{:});
        end
        
        %% Set Methods
        function this = set.X(this,definitionIn)
            this.X = this.setProperty(definitionIn);
        end
        
        function this = set.Y(this,definitionIn)
            this.Y = this.setProperty(definitionIn);
        end
        
        %% Parameterizable methods
        function this = refresh(this)
            warning('Point:refresh','Replace "refresh" with "rebuild."');
            this = rebuild(this);
        end
        
        function this = rebuild(this)
            this.X = rebuild(this.X);
            this.Y = rebuild(this.Y);
        end
    end
end