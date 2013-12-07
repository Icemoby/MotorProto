classdef RotatingMachineAssembly < CylindricalAssembly
    properties
        ElectricalFrequency = RotatingMachineAssembly.setProperty(60);
    end
    
    properties (Abstract,Dependent)
        AngularVelocity
    end
    
    properties (Dependent)
        GeometryMirrorSymmetry
    end
    
    methods
        function this = RotatingMachineAssembly(varargin)
            this = this@RotatingMachineAssembly(varargin{:});
        end

        function this = set.ElectricalFrequency(this,value)
            this.ElectricalFrequency = this.setProperty(value);
        end
        
        function boolOut = get.GeometryMirrorSymmetry(this)
            I = [this.Regions.IsUserDefined];
            if isempty(I) || all(~I)
                boolOut = true;
            else
                boolOut = this.Regions(I).hasPolarMirrorSymmetry;
            end
        end
        
      	function this = build(this,fraction)
            if nargin < 2
                if this.SolutionHalfWaveSymmetry
                    fraction = 1 / this.SolutionSpatialFrequency / 2;
                else
                    fraction = 1 / this.SolutionSpatialFrequency;
                end
            end
            
            %% Validate Inputs
            validateattributes(fraction,{'numeric'},{'scalar'});
            nCopies = this.GeometryFrequency * fraction;
            
          	assert(abs(nCopies-round(nCopies)) < sqrt(eps),...
                    'MotorProto:RotatingMachineAssembly:buildPitch',...
                    ['The input fraction must be an integer multiple'...
                    ' of the geometry frequencymust be an integer']);
              
            nCopies = round(nCopies);
            
            assert(nCopies > 0,...
                    'MotorProto:RotatingMachineAssembly:buildPitch',...
                    ['The GeometryFrequency times the input fraction'...
                    ' must be greater than one']);

            this = clean(this);
            
            regions = buildPreProcessing(this);
            
            %% Add copies of user input regions
            regions  = [regions,this.InputRegions];
            nRegions = numel(regions);
            if nRegions > 0
                regions                 = repmat(regions,1,nCopies);
                regions                 = copy(regions);
                [regions.IsUserDefined] = deal(false);
                
                angles = 2*pi / this.GeometryFrequency * ((1:nCopies) - 0.5);
                angles = repmat(angles,nRegions,1);
                angles = reshape(angles,nRegions*nCopies,1).';
                
                rotate(regions,angles,[0 0]);
                
                evenCopies = mod(0:(nCopies*nRegions-1), nRegions*2);
                evenCopies = (evenCopies > nRegions);
                
                if any(evenCopies)
                    reversePolarity(regions(evenCopies));
                end
                
                addModelRegion(this,regions);
            end
             
            %% Build the domain hull
            nRegions = nRegions * nCopies;
            ang      = 2 * pi * fraction;
            rot      = 0;
            rad      = [this.InnerRadius.Value,this.OuterRadius.Value];
          	dgOut    = Geometry2D.draw('Sector',...
                                        'Radius',rad,...
                                        'Angle',ang,...
                                        'Rotation',rot,...
                                        'PlotStyle',{'b'});
            this.DomainHull = copy(dgOut);
            if nRegions > 0
                dgOut = dgOut - [this.Regions.Geometry];
            end
            
            addModelRegion(this, [this.Name,'_DefaultRegion'],...
                                    dgOut,...
                                 	this.DefaultMaterial,...
                                   	'Static',...
                                   	-1);
            
            this = buildPostProcessing(this);
        end
                
        function regions = buildPreProcessing(this)
            %% Noop, subclasses may implement behavior here
            regions = [];
        end
        
        function this = buildPostProcessing(this)
            %% Noop, subclasses may implement behavior here
        end
        
      	function mesh = newMesh(this)
            mesh = RotatingMachineAssembly(this);
        end
    end
end