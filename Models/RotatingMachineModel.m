classdef RotatingMachineModel < Model
   	%RotatingMachineModel.m The standard class for objects representing Rotating Machines
    %
    % RotatingMachineModel methods:
    %   build - Populates the object's Assemblies property based on the current configuration
    %
	% RotatingMachineModel inheritance:
    %   RotatingMachineModel inherts methods and properties from Model.
    %
    %   See the help for Model for more information.
    %
    % See also Model, RotatingMachineMeshFactory, RotatingMachineAssembly
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$

    properties (Dependent, SetAccess = private)
        %% New Properties
        SpatialSymmetries
        
        %% Old Properties
        SolutionSpatialFrequency
        SolutionTemporalFrequency
        HasHalfWaveSymmetry
        SpaceModelFraction
    end
    
    methods
        %% Getters
        function value = get.SpatialSymmetries(this)
            assembly = this.Assemblies;
            value    = min([assembly.SpatialSymmetries]);
        end
        
        function value = get.SolutionSpatialFrequency(this)
            assemblies  = this.Assemblies;
            value       = assemblies(1).SpatialSymmetries;
            nAssemblies = numel(assemblies);
            for i = 2:nAssemblies
                value = gcd(value, assemblies(i).SpatialSymmetries);
            end
        end
        
        function value = get.SolutionTemporalFrequency(this)
            omega   = this.Assemblies.ElectricalFrequency;
            value   = [omega.Value];
            dF      = bsxfun(@minus, value.', value);
            isEqual = abs(dF) < sqrt(eps) * max(abs(value));
            if all(all(isEqual))
                value = value(1);
            end
        end
        
        function value = get.HasHalfWaveSymmetry(this)
            assemblies = this.Assemblies;
            value      = all([assemblies.HasHalfWaveSymmetry]);
        end
        
        function value = get.SpaceModelFraction(this)
            value = this.SolutionSpatialFrequency;
            if this.HasHalfWaveSymmetry
                value = value * 2;
            end
            value = 1 / value;
        end
        
        %% Build
        function this = build(this, symmetryType)
            symmetryType = 'space';
            % build - Populates the object's Assemblies property based on the current configuration
            % build(M) Constructs the geometry, mesh, and source information of
            % the model based on the current configuration of the objects in the
            % Assemblies array.  The build method attempts to construct the smallest
            % model possible based on a confluence of spatial and temporal
            % symmetries of all the entires in Assemblies property.
            %
            %   Example: Create a simple Synchronous IPM Machine
            %   S  = MotorProto('E.G.');
            %   M  = S.Model;
            %   %Create Rotor
            %   RT = M.newAssembly('myRotor', 'SynchronousRotor',...
            %                      'Poles', 4, 'InnerRadius', 0.25,...
            %                                  'OuterRadius', 0.5,...
            %                      'DefaultMaterial', IronExampleMaterial);
            %
            %  %Add Permanent Magnet
            %   PM = Geometry2D.draw('Rect', 'Width', 0.5, 'Length', 0.03,...
            %                        'Base', 'Center', 'Position', [0.35, 0]);
            %
            %   RT.addRegion('Magnet', PM,  PermanentMagnetExampleMaterial,...
            %                'Dynamic', 'Isolated');
            %
            %   %Create Stator
            %   ST = M.newAssembly('myStator', 'Stator',...
            %                      'Poles', 4, 'Teeth', 24, 'InnerRadius', 0.51,...
            %                      'OuterRadius', 1,...
            %                      'DefaultMaterial', IronExampleMaterial);
            %
            %   %Add a slot using template
            %   [slotBody, slotFront] = slotTemplate(24,  0.51, 1, 0.05, 0.01, ...
            %                                        0.4, 0.5,  1, 'Auto');
            %
            %   ST.Slot.Shape = slotBody;
            %   ST.Slot.Turns = 2;
            %
            %   ST.addRegion('SlotFront', slotFront, Air, 'Static', 'Isolated');
            %
            %   figure;
            %   title('Preview');
            %   M.preview;
            %
            %   M.build;
            %   figure;
            %   title('One Pole Plot After Calling build(M)');
            %   M.plot;
            %
            %   figure;
            %   title('Mesh with default settings');
            %   M.Mesh.plot;
            %
            % See also MotorProto, RotatingMachineModel

            if nargin == 1 || isempty(symmetryType)
                symmetryType = 'space';
            end
            
            assemblies  = this.Assemblies;
          	nAssemblies = numel(assemblies);
            switch lower(symmetryType)
                case 'space'
                    fraction = 1 / this.SpatialSymmetries;
                    if this.HasHalfWaveSymmetry
                        fraction = fraction / 2;
                    end
                    for i = 1:nAssemblies
                        assemblies(i).build(fraction);
                    end
                case 'spacetime'
                    warning('MotorProto:Verbose', 'This method is not completely general');
                    for i = 1:nAssemblies
                        assemblies(i).build(1 / assemblies(i).SolutionSpaceTimeSymmetry(1));
                    end
                otherwise
                    error('MotorProto:RotatingMachineModel', 'Unknown symmetryType %s', lower(symmetryType));
            end
            
            warning('MotorProto:Verbose', 'Separate out this call');
            build(this.Mesh);
        end
    end
end