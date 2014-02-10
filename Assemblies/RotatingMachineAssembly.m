classdef RotatingMachineAssembly < CylindricalAssembly
    %RotatingMachineAssembly.m An abstract class representing Rotating Machines assemblies
    %   Subclasses of the RotatingMachineAssembly represent high level logical
    %   units of electric machines such as stators and rotors.
    %
    % RotatingMachineAssembly methods:
    %   build - Constructs the object's Components array based on the current configuration
    %
    % RotatingMachineAssembly properties:
    %   ElectricalFrequency - The fundamental electrical frequency of the object
    %   InitialAngle        - The angle of the RotatingMachineAssembly object at time zero
    %
    % RotatingMachineAssembly inherits properties and methods CylindricalAssembly.
    %
    % See the help files for CylindricalAssembly for more information.
    %
    % See also MotorProto, Model, CylindricalAssembly, Stator, SynchronousRotor
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%ElectricalFrequency - The fundamental electrical frequency of the object
    %   The ElectricalFrequency property specifies the fundamental electrical
    %   frequency of the RotatingMachineAssembly object.
    %
    % See also RotatingMachineAssembly
    ElectricalFrequency;
    
 	%InitialAngle - The angle of the RotatingMachineAssembly object at time zero
    %   The InitialAngle property specifies the angle (in radians) that the
    %   RotatingMachineAssembly object is rotated in the counter-clockwise
    %   direction at time zero when measured from the x-axis.
    %
    % See also RotatingMachineAssembly
    InitialAngle;
%}

    properties
        ElectricalFrequency = RotatingMachineAssembly.setProperty(60);
        InitialAngle        = RotatingMachineAssembly.setProperty(0);
    end
    
    properties (Abstract, Dependent)
        AngularVelocity
    end
    
    properties (SetAccess = protected)
        ModeledFraction
    end
    
    properties (Dependent)
        GeometryMirrorSymmetry
    end
    
    methods
        %% Constructor
        function this = RotatingMachineAssembly(varargin)
            this = this@CylindricalAssembly(varargin{:});
        end

        %% Setters
        function this = set.ElectricalFrequency(this, value)
            this.ElectricalFrequency = this.setProperty(value);
        end
        
        function this = set.InitialAngle(this, value)
            this.InitialAngle = this.setProperty(value);
        end
        
        %% Getters
        function boolOut = get.GeometryMirrorSymmetry(this)
            boolOut = this.InputRegions.hasPolarMirrorSymmetry;
        end
        
        %% Others
      	function this = build(this, modeledFraction)
            %build - Populates the object's Components property based on the current configuration
            %build(A) uses the user defined regions in the InputRegions property
            %along with other object properties to populate the Components
            %property array. The build method attempts to construct the smallest
            %model possible based on the number of spatial symmetries.
            %
            %   Example: Create a simple Stator
            %   S = MotorProto('E.G.');
            %   M = S.Model;
            %
            %   %Create Stator
            %   A = M.newAssembly('myStator', 'Stator',...
            %                      'Poles', 4, 'Teeth', 24, 'InnerRadius', 0.51,...
            %                      'OuterRadius', 1,...
            %                      'DefaultMaterial', IronExampleMaterial);
            %
            %   %Add a slot using template
            %   [slotBody, slotFront] = slotTemplate(24,  0.51, 1, 0.05, 0.01, ...
            %                                        0.4, 0.5,  1, 'Auto');
            %
            %   A.Slot.Shape = slotBody;
            %   A.Slot.Turns = 2;
            %
            %   A.addRegion('SlotFront', slotFront, Air, 'Static', 'Isolated');
            %
            %   figure;
            %   title('Preview');
            %   A.preview;
            %
            %   A.build;
            %   figure;
            %   title('One Pole Plot After Calling build(A)');
            %   A.plot;
            %   
            %   figure;
            %   title('Wireframe');
            %   A.wireframe;
            %
            % See also MotorProto, Model, RotatingMachineAssembly
            
            if nargin < 2
                if this.HasHalfWaveSymmetry
                    modeledFraction = 1 / this.SpatialSymmetries / 2;
                else
                    modeledFraction = 1 / this.SpatialSymmetries;
                end
            end
            
            %% Validate Inputs
            validateattributes(modeledFraction, {'numeric'}, {'scalar'});
            nCopies = this.GeometricSymmetries * modeledFraction;
            
          	assert(abs(nCopies-round(nCopies)) < sqrt(eps), 'MotorProto:RotatingMachineAssembly', 'ModelFraction must be an integer multiple of 1 / GeometricSymmetries');
              
            nCopies = round(nCopies);
            
            assert(nCopies > 0, 'MotorProto:RotatingMachineAssembly', 'GeometricSymmetries * ModeledFraction >= 1 is required');

            this = clean(this);
            
            [conductors, nonConductors, connectionMatrix] = buildPreProcessing(this);
            
            %% Add copies of user input regions
            regions  = [conductors, nonConductors, this.InputRegions];
            nRegions = numel(regions);
            if nRegions > 0
                %% Copy the regions and rotate them
                regions                 = repmat(regions, 1, nCopies);
                regions                 = copy(regions);
                [regions.IsUserDefined] = deal(false);
                
                angles = 2 * pi / this.GeometricSymmetries * ((1:nCopies) - 0.5) + this.InitialAngle.Value;
                angles = repmat(angles, nRegions,1);
                angles = reshape(angles, nRegions * nCopies, 1).';
                
                rotate(regions, angles, [0 0]);
                
                evenCopies = mod(0:(nCopies * nRegions-1), nRegions * 2)+1;
                evenCopies = (evenCopies > nRegions);
                
                if any(evenCopies)
                    reversePolarity(regions(evenCopies));
                end
                
                addModelRegion(this, regions);
                
                %% Copy the connection matrix and renumber
                if ~isempty(connectionMatrix)
                    offset     = ((0:(nCopies-1)) * nRegions);
                    newMatrix  = connectionMatrix;
                    for i = 2:nCopies
                        newMatrix = newMatrix + (offset(i) - offset(i-1)) * sign(newMatrix);
                        connectionMatrix = cat(3,connectionMatrix,newMatrix);
                    end
                end
                
                %% Determine how the regions were stored in the assembly
                %   Loop renumbers the connection matrix based on the
                %   mapping between the conductors/regions and their
                %   location in the components array after running the
                %   addModelRegion command
                
                components = this.Components;
                for i = 1:size(connectionMatrix,1)
                    for j = 1:size(connectionMatrix,2)
                        for k = 1:size(connectionMatrix,3)
                            l = find(regions(connectionMatrix(i,j,k)) == components);
                            connectionMatrix(i,j,k) = l;
                        end
                    end
                end
            end
             
            %% Build the domain hull
            nRegions = nRegions * nCopies;
            ang      = 2 * pi * modeledFraction;
            rot      = this.InitialAngle.Value;
            rad      = [this.InnerRadius.Value, this.OuterRadius.Value];
          	dgOut    = Geometry2D.draw('Sector', 'Radius', rad, 'Angle', ang, 'Rotation', rot, 'PlotStyle', {'b'});
            this.DomainHull = copy(dgOut);
            if nRegions > 0
                dgOut = dgOut - [this.Regions.Geometry];
            end
            
            addModelRegion(this, [this.Name,'_DefaultRegion'], dgOut, this.DefaultMaterial, DynamicsTypes.Static, -1);
            
            this                 = buildPostProcessing(this, connectionMatrix, modeledFraction);
            this.ModeledFraction = modeledFraction;
        end
                
        function [conductors, nonConductors, connectionMatrix] = buildPreProcessing(this)
            %% Noop, subclasses may implement behavior here
            conductors       = [];
            nonConductors    = [];
            connectionMatrix = [];
        end
        
        function this = buildPostProcessing(this, connectionMatrix, modeledFraction)
            %% Noop, subclasses may implement behavior here
        end
        
      	function mesh = newMesh(this)
            mesh = RotatingMachineMeshFactory(this);
        end
    end
    
    methods (Sealed)        
        function [t, h] = getTimePoints(this, Nt)
            nAssemblies = numel(this);
            
            %% Calculate Electrical Frequency
            f           = [this.ElectricalFrequency];
            f           = [f.Value];
            assert(all(abs(f - f(1))) < sqrt(eps) * f(1), 'MotorProto:RotatingMachineAssembly', 'All assembly electrical frequencies must be the same');
            T           = 1 / mean(f);
            
            %% Calculate Fundamental Period (possibly subharmonic)
            Np = zeros(1,nAssemblies);
            for i = 1:nAssemblies
                Np(i) = this(i).Poles.Value;
            end
            assert(all(abs(Np - Np(1)) < sqrt(eps) * Np(1)));
            Np = Np(1);
            
            Mf = max([this.ModeledFraction]);
            
            if Np * Mf > 2
                T  = T  * ceil(Np * Mf / 2);
                Nt = Nt * ceil(Np * Mf / 2); 
            end
            
            Nh          = ceil(Nt / 2);
            NhSynch     = floor((Nh + 2) / 6) * 6;
            NhStat      = NhSynch + 1;
            hStat       = 1:2:NhStat;
            hSynch      = 0:6:NhSynch;
            Nt          = (2 * NhSynch + 6) + 1;
            t           = linspace(0, 1, Nt) * T;
            h           = cell(1, nAssemblies);
            
            for i = 1:nAssemblies
                if this(i).AngularVelocity == 0
                    h{i} = hStat;
                else
                    h{i} = hSynch;
                end
            end
        end
    end
end