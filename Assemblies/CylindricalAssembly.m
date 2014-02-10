classdef CylindricalAssembly < Assembly
    %CylindricalAssembly.m An abstract class representing Assembly objects defined on a cylindrical annulus
    %   Subclasses of the CylindricalAssembly class represent Assembly objects
    %   which have a DomainHull defined on a cylindrical annulus.
    %
    % CylindricalAssembly properties:
    %   InnerRadius - The inner radius of the CylindricalAssembly
    %   OuterRadius - The outer radius of the CylindricalAssembly
    %   Length      - The length of the CylindricalAssembly
    %
    % CylindricalAssembly inherits properties and methods from Assembly.
    %
    % See the help file for Assembly for more information.
    %
    % See also MotorProto, Model, Assembly, 
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%InnerRadius - The inner radius of the CylindricalAssembly object
    %
    % See also CylindricalAssembly
    InnerRadius;
    
 	%OuterRadius - The outer radius of the CylindricalAssembly object
    %
    % See also CylindricalAssembly
    OuterRadius;
    
 	%Length - The length of the CylindricalAssembly object
    %
    % See also CylindricalAssembly
    Length;
%}
    properties
        %% Basic Geometry Properties
        InnerRadius = CylindricalAssembly.setProperty(0)
        OuterRadius = CylindricalAssembly.setProperty(1)
        Length      = CylindricalAssembly.setProperty(1)
    end
    
	properties (Abstract,Dependent)
        %% New Symmetry Properties
        SpatialSymmetries
        HasHalfWaveSymmetry
        SpaceTimeSymmetries

        GeometricSymmetries
        
        %% Old Symmetry Properties
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        SolutionSpaceTimeSymmetry
        SolutionSpaceTimeCoefficients
        GeometryFrequency
        GeometryMirrorSymmetry
    end
    
    properties (Dependent, SetAccess = private)
        Mass
        TangentialBoundaries
    end
    
    methods
        %% Constructions
     	function this = CylindricalAssembly(varargin)
            this = this@Assembly(varargin{:});
        end
        
        %% Setters
        function set.InnerRadius(this,radiusIn)
            this.InnerRadius = this.setProperty(radiusIn);
        end
        
        function set.OuterRadius(this,radiusIn)
            this.OuterRadius = this.setProperty(radiusIn);
        end
        
        function set.Length(this,lengthIn)
            this.Length = this.setProperty(lengthIn);
        end
        
        %% Getters
        function angles = get.TangentialBoundaries(this)
            if ~isempty(this.DomainHull)
                da     = this.DomainHull.Angle.Value;
                a      = this.DomainHull.Rotation;
                angles = [a,a + da];
            else
                angles = [];
            end
        end
        
        function mass = get.Mass(this)
            regions  = this.Regions;
            nRegions = length(regions);
            mass     = 0;
            for i = 1:nRegions
                mass = mass + regions(i).Geometry.area * regions(i).Material.Density;
            end
            mass = mass * this.Length.Value;
            mass = mass / this.ModeledFraction;
        end
        
        %% Element Functions
        function previewElement(this)
            dh = makeElementDomainHull(this);
            
            if ~isempty(this.InputRegions);
                inputGeometry = [this.InputRegions.Geometry];
                dh = dh - inputGeometry;
            else
                inputGeometry = [];
            end
            
            plot([dh, inputGeometry]);
        end
        
        function hull = makeElementDomainHull(this)
            r    = [this.InnerRadius.Value, this.OuterRadius.Value];
            a    = 2 * pi / this.GeometricSymmetries;
            hull = Geometry2D.draw('Sector', 'Radius', r, 'Angle', a, 'Rotation', -a/2);
        end
    end
end