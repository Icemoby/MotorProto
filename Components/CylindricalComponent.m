classdef CylindricalAssembly < Assembly
    properties
        %% Basic Geometry Properties
        InnerRadius = CylindricalAssembly.setProperty(0)
        OuterRadius = CylindricalAssembly.setProperty(1)
        Length      = CylindricalAssembly.setProperty(1)
        
        %% Material Properties
        DefaultMaterial
    end
    
	properties (Abstract,Dependent)
        SolutionSpatialFrequency
        SolutionHalfWaveSymmetry
        SolutionTemporalFrequency
        SolutionTemporalSymmetry
        SolutionSpaceTimeSymmetry
        GeometryFrequency
        GeometryMirrorSymmetry
    end
    
    properties (Dependent,SetAccess = private)
        TangentialBoundaries
    end
    
    methods
        function angles = get.TangentialBoundaries(this)
            if ~isempty(this.DomainHull)
                da     = this.DomainHull.Angle.Value;
                a      = this.DomainHull.Rotation;
                angles = [a,a + da];
            else
                angles = [];
            end
        end
        
     	function this = CylindricalAssembly(varargin)
            this = this@Assembly(varargin{:});
        end
        
        function set.InnerRadius(this,radiusIn)
            this.InnerRadius = this.setProperty(radiusIn);
        end
        
        function set.OuterRadius(this,radiusIn)
            this.OuterRadius = this.setProperty(radiusIn);
        end
        
        function set.Length(this,lengthIn)
            this.Length = this.setProperty(lengthIn);
        end
        
       	function set.DefaultMaterial(this,materialIn)
            assert(isa(materialIn,'MaterialProperty'),...
                    'MotorProto:CylindricalComponent:invalidType',...
                    'DefaultMaterial must be a MaterialProperty object');
            this.DefaultMaterial = materialIn;
        end
        
        function previewElement(this)
            dh = makeElementDomainHull(this);
            dh = dh - [this.InputRegions.Geometry];
            plot([dh,this.InputRegions.Geometry]);
        end
        
        function hull = makeElementDomainHull(this)
            r    = [this.InnerRadius.Value,this.OuterRadius.Value];
            a    = 2 * pi / this.GeometryFrequency;
            hull = Geometry2D.draw('Sector','Radius',r,....
                                            'Angle',a,....
                                            'Rotation',-a/2);
        end
    end
end