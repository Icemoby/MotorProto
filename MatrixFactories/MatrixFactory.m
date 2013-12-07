classdef MatrixFactory
    properties (SetAccess = protected,Dependent)
        %% Handle objects derived from Model_, copy on return
        Model
        Mesh
        Assemblies
    end
    
    properties (SetAccess = protected)
        %% Matrix components are stored in structures for fast access
        
        %% Discrete del operator matrices, div, grad, curl
        Del
        
        %% Boundary matrices
        Boundary
        
        %% M(t)*x_t + K(t)*x + g(x,t) = F*f(t);
        Mass
        Stiffness
        Jacobian
        Exogenous

        %% Post Processing
        PostProcessing    
        
        %% Indices
        Index
    end
    
    properties (Hidden,SetAccess = protected)
        %% Input Argument
        Model_
    end
    
    properties (Hidden,SetAccess = protected,Dependent)
        %% Internal arguments, no copy on output
        Mesh_
        Assemblies_
    end
    
    methods
        %% Move These To StaticMatrixFactory
        %% Constructor
        function this = MatrixFactory(model)
            warning('MotorProto:Verbose',...
                    'Present implementation assumes anti-periodic symmetry');
            if nargin > 0
                this.Model  = model;
            end
        end

        %% Setters
        function this = set.Model(this,model)
            this.Model_  = model;
        end
        
        %% Getters
        function value = get.Model(this)
            value = copy(this.Model_);
        end
        
        function value = get.Mesh(this)
            value = copy(this.Mesh_);
        end
        
        function value = get.Mesh_(this)
            value = this.Model_.Mesh;
        end
        
        function value = get.Assemblies(this)
            value = copy(this.Assemblies_);
        end
        
        function value = get.Assemblies_(this)
            value = [this.Mesh_.Assembly];
        end
        
        %% preprocessing functions
        function this  = build(this, symmetryType)
            build(this.Model_, symmetryType);

            this = buildIndexVectors(this);
            
            this = buildDelMatrices(this);
            
            this = buildPostProcessingMatrices(this);
            
            this = buildBoundaryMatrices(this);
            
            this = buildMassMatrices(this);
            this = buildStiffnessMatrices(this);
            this = buildJacobianMatrices(this);
            this = buildExogenousMatrices(this);
            
            this = applyBoundaryConditions(this);
        end
        
        function this  = buildIndexVectors(this)
            this.Index = buildMatrices(this, this.Index, 'buildLocalIndexVectors');
            this.Index = buildGlobalIndexVectors(this, this.Index);
        end
        
        function index = buildGlobalIndexVectors(this,index)
            index.Global(1)   = index.Local(1);
            index.Local(1).X  = 1:index.Local(1).Unknowns;
            index.Global(1).X = 1:index.Global(1).Unknowns;
            nUnknowns         = index.Local(1).Unknowns;
            
            for i = 2:numel(index.Local);
             	index.Global(i)   = buildGlobalIndexField(this, index.Local(i), nUnknowns);
                index.Local(i).X  = 1:index.Local(i).Unknowns;
                index.Global(i).X = (nUnknowns + 1):index.Global(i).Unknowns;
                nUnknowns         = index.Global(i).Unknowns;
            end
        end
        
        function structure = buildGlobalIndexField(this, structure, nUnknowns)
            nStructs   = numel(structure);
            fieldNames = fields(structure);
            
            for i = 1:nStructs
                for j = 1:numel(fieldNames)
                    if isnumeric(structure(i).(fieldNames{j}))
                        structure(i).(fieldNames{j}) = structure(i).(fieldNames{j}) + nUnknowns;
                    elseif iscell(structure(i).(fieldNames{j}))
                        n = numel(structure(i).(fieldNames{j}));
                        for k = 1:n
                            structure(i).(fieldNames{j}){k} = structure(i).(fieldNames{j}){k} + nUnknowns;
                        end
                    else
                        structure(i).(fieldNames{j})  = buildGlobalIndexField(this, structure(i).(fieldNames{j}), nUnknowns);
                    end
                end
            end
        end
        
        function [structure, newFields, pbcType] = buildLocalIndexVectors(this, structure, iMesh)
            newFields = {'Local'};
            pbcType   = {'None'};
            
            if isfield(structure,'Local')
                local = structure.Local;
            else
                local = struct([]);
            end
            
            mesh   = this.Mesh_(iMesh);
            nNodes = numel(mesh.X);
            
            local(iMesh).A        = 1:nNodes;
            local(iMesh).Unknowns = nNodes;
            
            local(iMesh).Boundary.Tangent(1).Nodes = mesh.PeriodicBoundaryNodes(1,:);
            local(iMesh).Boundary.Tangent(2).Nodes = mesh.PeriodicBoundaryNodes(1,:);
            
            local(iMesh).Boundary.Radius(1).Nodes  = unique(mesh.RadialBoundaryEdges{1});
            local(iMesh).Boundary.Radius(2).Nodes  = unique(mesh.RadialBoundaryEdges{2});
            
            structure.Local = local;
        end
        
        %% Del Matrices
        function this = buildDelMatrices(this)
            this.Del = buildMatrices(this, this.Del, 'buildCurlMatrices');
        end
        
        function [structure, newFields, pbcType] = buildCurlMatrices(this, structure, iMesh)
            %% New Field Definitions
            newFields = {'Curl'};
            pbcType   = {'None'};
            
            %% Build Matrices
            mesh      = this.Mesh_(iMesh);
           	el        = mesh.Elements;
            elArea    = mesh.ElementAreas;
            nElements = length(el);
            nUnknowns = this.Index.Local(iMesh).Unknowns;
            
        	x = mesh.X(el);
            y = mesh.Y(el);
            b =  [y(2,:)-y(3,:); y(3,:)-y(1,:); y(1,:)-y(2,:)];
            c = -[x(2,:)-x(3,:); x(3,:)-x(1,:); x(1,:)-x(2,:)];
            
            i = [1:nElements,1:nElements,1:nElements];
            j = [el(1,:),el(2,:),el(3,:)];
            
            %% Make node to element curl matrices
            s = -[b(1,:)./elArea, b(2,:)./elArea, b(3,:)./elArea] / 2;
            structure.Curl(iMesh).Zn2Ye = sparse(i,j,s,nElements,nUnknowns);
            
            s = [c(1,:)./elArea, c(2,:)./elArea, c(3,:)./elArea] / 2;
            structure.Curl(iMesh).Zn2Xe = sparse(i,j,s,nElements,nUnknowns);

            %% Make element to node curl matrices
            s = [b(1,:), b(2,:), b(3,:)] / 2;
            structure.Curl(iMesh).Ye2Zn = sparse(j,i,s,nUnknowns,nElements);

            s = -[c(1,:), c(2,:), c(3,:)] / 2;
            structure.Curl(iMesh).Xe2Zn = sparse(j,i,s,nUnknowns,nElements);
        end
        
        %% Boundary Matrices
        function this = buildBoundaryMatrices(this)
            this.Boundary.Radial = buildMatrices(this, struct(''), 'buildRadialBoundaryMatrices');
        end
        
        function [structure, newFields, pbcType] = buildRadialBoundaryMatrices(this, structure, iMesh)
            %% New Field Definitions
            newFields = {'R'   , 'D'                , 'S'  , 'P'  , 'F'   ,'G'};
            pbcType   = {'None', 'BoundaryOperation', 'Row', 'Row', 'None', 'None'};
            
            if isempty(structure)
                structure = cell2struct(cell(1, numel(newFields)), newFields, 2);
            end
            
            %% Get boundary radii
            meshes    = this.Mesh_;
            nMeshes   = numel(meshes);
            r         = [meshes.RadialBoundaryRadii];
            rSelf     = r((2*iMesh-1):(2*iMesh));
            rOpposite = [0 inf];
            
            if iMesh > 1
                rOpposite(1) = r(2*iMesh - 2);
            end
            
            if iMesh < nMeshes
                rOpposite(2) = r(2*iMesh + 1);
            end
            
            %% get boundary edges
            mesh       = meshes(iMesh);
            bEdges     = mesh.RadialBoundaryEdges;
            
            %% get angular velocity
            assembly   = mesh.Assembly;
            bVelocity  = assembly.AngularVelocity;
            bHarmonics = this.getRadialBoundaryHarmonics(iMesh);
            
            %% determine modeled fraction
            mFraction  = assembly.ModeledFraction;
            
            %% calculate matrices
            nUnknowns  = this.Index.Local(iMesh).Unknowns;
            for i = 1:2
                n          = bHarmonics{i}.';
                nHarmonics = numel(n);
            
                if rOpposite(i) == 0 || isinf(rOpposite(i))
                    R = sparse(1, 1, 1);
                else
                    R = bVelocity(1) * n;
                end
            
                x = mesh.X(bEdges{i});
                y = mesh.Y(bEdges{i});
                t = atan2(y, x);

                [t,K] = sort(t);
                b     = bEdges{i};
                for k = 1:length(K)
                    b(:, k) = b(K(:,k), k);
                end
            
                %% dft matrix
                s1 =   1i * bsxfun(@rdivide, exp(-1i * n * t(2,:)), n) + (exp(-1i * n * t(1,:)) - exp(-1i * n * t(2,:))) ./ (n.^2 * (t(1,:) - t(2,:)));

                s2 =  -1i * bsxfun(@rdivide, exp(-1i * n * t(1,:)), n) + (exp(-1i * n * t(2,:)) - exp(-1i * n * t(1,:))) ./ (n.^2 * (t(1,:) - t(2,:)));

%                 s  = [ reshape(s1, 1, []), reshape(s2, 1, [])] / mFraction / 2 / pi;
                s  = [ reshape(s1, 1, []), reshape(s2, 1, [])] / 2 / pi;            
                k  = [ reshape((1:nHarmonics).' * ones(size(b(2, :))),1,[]), reshape((1:nHarmonics).' * ones(size(b(1,:))), 1, [])];
                    
                l  = [ reshape(ones(size(n))*b(2,:), 1, []), reshape(ones(size(n))*b(1,:), 1, [])];
            
            	D = sparse(k,l,s,nHarmonics,nUnknowns);
%                 S = (-1)^(i) * D' * 2 * pi * rSelf(i) * mFraction;
                S = (-1)^(i) * D' * 2 * pi * rSelf(i);
                %% idft matrix
                [k,getUnique] = unique(b);
                t             = t(getUnique).';

                k = repmat(k, 1, nHarmonics);
                l = repmat(1:nHarmonics, numel(t), 1);
                s = exp(1i * (n * t).');

                k = reshape(k,1,[]);
                l = reshape(l,1,[]);
                s = reshape(s,1,[]);
            
                P = sparse(k, l, s, nUnknowns,nHarmonics);
            
                %% boundary transfer matrices
                if isinf(rOpposite(i))
                    sf =   abs(n) / rSelf(i);
                    sg =   zeros(1,nHarmonics);
                elseif rOpposite(i) == 0;
                    sf = - abs(n) / rSelf(i);
                    sg =   zeros(1,nHarmonics);
                else
                    a  = rOpposite(i) / rSelf(i);
                    sf = (n     .* (a.^n + a.^(-n))) ./ (rSelf(i)     * (a.^n - a.^(-n)));
                    sg = (2 * n                    ) ./ (rOpposite(i) * (a.^n - a.^(-n)));
                end
                F = sparse(1:nHarmonics, 1:nHarmonics, sf / mu_o, nHarmonics, nHarmonics);
                G = sparse(1:nHarmonics, 1:nHarmonics, sg / mu_o, nHarmonics, nHarmonics);
                
                if i == 1
                    structure.R(iMesh).Inner = R;
                    structure.D(iMesh).Inner = D;
                    structure.S(iMesh).Inner = S;
                    structure.P(iMesh).Inner = P;
                    structure.F(iMesh).Inner = F;
                    structure.G(iMesh).Inner = G;
                else
                    structure.R(iMesh).Outer = R;
                    structure.D(iMesh).Outer = D;
                    structure.S(iMesh).Outer = S;
                    structure.P(iMesh).Outer = P;
                    structure.F(iMesh).Outer = F;
                    structure.G(iMesh).Outer = G;
                end
            end
        end
        
        %% Stiffness Matrices
        function this = buildStiffnessMatrices(this)
        	this.Stiffness = buildMatrices(this, this.Stiffness, 'buildReluctivityMatrices');
        end
        
      	function [structure, newFields, pbcType] = buildReluctivityMatrices(this, structure, iMesh)
            %% New Field Definitions
            newFields = {'Reluctivity'};
            pbcType   = {'Both'};
            
            if isfield(structure, newFields{1})
                reluctivity = structure.Reluctivity;
            else
                reluctivity = struct([]);
            end
            
            %% Build Matrices
            curl                     = this.Del.Curl(iMesh);
            reluctivity(iMesh).El2El = curl.Xe2Zn*curl.Zn2Xe + curl.Ye2Zn*curl.Zn2Ye;
            reluctivity(iMesh).El2El = - reluctivity(iMesh).El2El / mu_o;
            
            %% Append to Structure
            structure.Reluctivity = reluctivity;
        end
        
        %% Mass Matrices
        function this = buildMassMatrices(this)
            this.Mass = buildMatrices(this, this.Mass, 'buildConductivityMatrices');
        end
        
        function [structure, newFields, pbcType] = buildConductivityMatrices(this, structure, iMesh)
            %% New Field Definitions
          	newFields = {'Conductivity'};
            pbcType  = {'Both'};
            
            if isfield(structure, newFields{1})
                conductivity = structure.Conductivity;
            else
                conductivity = struct([]);
            end
            
            %% Build Matrices
            mesh      = this.Mesh_(iMesh);
            regions   = mesh.Regions;
            materials = [regions.Material];
            elCond    = [materials.sigma];
            
            isStatic  = regions.hasStaticDynamics;
            elRegions = mesh.ElementRegions;
            I         = ismember(elRegions, find(~isStatic));
            
            elRegions = elRegions(I);
            el        = mesh.Elements(:,I);
            elArea    = mesh.ElementAreas(I);
            elCond    = elCond(elRegions);
            nUnknowns = this.Index.Local(iMesh).Unknowns;
            
            %% Make element to element coupling (through A_t)
            i = [el(1,:), el(2,:), el(3,:),...
                 el(1,:), el(2,:), el(1,:),...
                 el(2,:), el(3,:), el(3,:)];

            j = [el(1,:), el(2,:), el(3,:),...
                 el(2,:), el(3,:), el(3,:),...
                 el(1,:), el(2,:), el(1,:)];
             
            s = elCond .* elArea / 12;
            s = [repmat(2*s,1,3), repmat(s,1,6)];
            
            conductivity(iMesh).El2El = sparse(i,j,s,nUnknowns,nUnknowns);
            
            %% Append to Structure
            structure.Conductivity = conductivity;
        end
        
        %% Jacobian Matrices
        function this = buildJacobianMatrices(this)
            this.Jacobian = buildMatrices(this, this.Jacobian, 'buildMagnetizationJacobianMatrices');
            this.Jacobian = buildMatrices(this, this.Jacobian, 'buildFluxDensityJacobianMatrices'  );
        end
        
        function [structure, newFields, pbcType] = buildMagnetizationJacobianMatrices(this, structure, iMesh)
            %% New Field Definitions
            newFields = {'MagnetizationCurrent'};
            pbcType  = {'Row'};
            
            if isfield(structure, newFields{1})
                magnetizationCurrent = structure.MagnetizationCurrent;
            else
                magnetizationCurrent = struct([]);
            end
            
            %% Build Matrices
            curl = this.Del.Curl(iMesh);
            
            magnetizationCurrent(iMesh).dIzdMx = curl.Xe2Zn;
            magnetizationCurrent(iMesh).dIzdMy = curl.Ye2Zn;
            
            %% Update Structure
            structure.MagnetizationCurrent = magnetizationCurrent;
        end
        
        function [structure, newFields, pbcType] = buildFluxDensityJacobianMatrices(this, structure, iMesh)
            %% Neew Field Definitions
            newFields = {'FluxDensity'};
            pbcType   = {'Column'};
            
            if isfield(structure, newFields{1})
                fluxDensity = structure.FluxDensity;
            else
                fluxDensity = struct([]);
            end
            
            %% Build Matrices
            curl = this.Del.Curl(iMesh);
            
            fluxDensity(iMesh).dBxdXz = curl.Zn2Xe;
            fluxDensity(iMesh).dBydXz = curl.Zn2Ye;
            
            %% Update Structure
            structure.FluxDensity = fluxDensity;
        end
        
        %% Exogenous Input Matrices
        function this = buildExogenousMatrices(this)
            this.Exogenous = buildMatrices(this, this.Exogenous, 'buildMagneticInputMatrices');
        end
        
        function [structure, newFields, pbcType] = buildMagneticInputMatrices(this, structure, iMesh)
            %% New Field Definitions
            newFields = {'Magnetic'};
            pbcType   = {'Row'};
            
            if isfield(structure, newFields{1})
                mqsInput = structure.Magnetic;
            else
                mqsInput = struct([]);
            end
            
            warning('MotorProto:Verbose', 'Put the conductivity/length into the source/region coefficients');
            warning('MotorProto:Verbose', 'This function might be made faster. E.g. short circuit for no sources');
            %error('Update for new source class, multiple parallel paths, etc. Calculate total current from each parallel path as input to take into account unequal resistances due to unequal geometry discretization');
            
            %% Build Matrices
            row = zeros(1,0);
            col = zeros(1,0);
            val = zeros(1,0);
            
            mesh       = this.Mesh_(iMesh);
            assembly   = mesh.Assembly;
            regions    = assembly.Regions;
            sources    = assembly.Sources;
            components = assembly.Components;
            geometry   = [regions.Geometry];
            
            nComponents = numel(components); 
            nRegions    = numel(regions);
            nSources    = numel(sources);
            nUnknowns   = this.Index.Local(iMesh).Unknowns;
            
            %% Create Current Density to Nodal Current Conversion Matrix
            if nSources > 0
                ind = 0;
                nEntries = numel(mesh.Elements);
                rows     = zeros(1, nEntries);
                cols     = zeros(1, nEntries);
                vals     = zeros(1, nEntries);
                for i = 1:nComponents
                    jRegion = find(components(i) == regions);
                    if ~isempty(jRegion)
                        iElement = (mesh.ElementRegions == jRegion);
                        iNode    = reshape(       mesh.Elements(:,  iElement)       , 1, []);
                        vAreas   = reshape(repmat(mesh.ElementAreas(iElement), 3, 1), 1, []) / 3;
                        iElement = reshape(repmat(             find(iElement), 3, 1), 1, []);
                        
                        I        = (ind+1):(ind+numel(iNode));
                        
                        rows(I)  = iNode;
                        cols(I)  = iElement;
                        vals(I)  = vAreas;
                        
                        ind      = ind + numel(iNode);
                    end
                end
                J2Node = sparse(rows, cols, vals, nUnknowns, length(mesh.Elements));
            end

%             %% Create Source to Current Density Conversion Matrix
%             for i = 1:nSources
%                 nPhases            = sources(i).Phases.Value;
%                 connectionMatrices = sources(i).ConnectionMatrices;
%                 connectionPolarity = sources(i).ConnectionPolarity;
%                 
%                 phaseResistance    = zeros(1,nPhases);
%                 maxEntries         = 0;
%                 
%                 %% Calculate Resistances
%                 for j = 1:nPhases
%                     [nSeries, nParallel]   = size(connectionMatrices{j});
%                     maxEntries               = maxEntries + nSeries * nParallel * nPhases;
%                     connectionRegions      = components(connectionMatrices{j});
%                     connectionGeometry     = [connectionRegions.Geometry];
%                     connectionMaterial     = [connectionRegions.Material];
%                     connectionConductivity = [connectionMaterial.Conductivity];
%                     connectionArea         = [connectionGeometry.area];
%                     connectionResistance   = assembly.Length.Value ./ (connectionConductivity .* connectionArea);
%                     connectionResistance   = reshape(connectionResistance, nSeries, nParallel);
%                     
%                     seriesResistance       = sum(connectionResistance, 1);
%                     phaseResistance(j)     = 1 ./ sum(1 ./ seriesResistance, 2);
%                 end
%                 parallelResistance = 1 ./ sum(1 ./ phaseResistance, 2);
%                 
%                 %% Calculate Matrix
%                 if (sources(i).ConnectionType == ConnectionTypes.Wye)
%                     rows = zeros(1, maxEntries);
%                     cols = zeros(1, maxEntries);
%                     vals = zeros(1, maxEntries);
%                     
%                     for j = 1:nPhases
%                         ind  = 0;
%                         [nSeries, nParallel]   = size(connectionMatrices{j});
%                         connectionRegions      = components(connectionMatrices{j});
%                         connectionMaterial     = [connectionRegions.Material];
%                         connectionConductivity = [connectionMaterial.Conductivity];
%                         
%                         if (sources(i).Type == SourceTypes.VoltageSource) && (assembly.ConnectionType == ConnectionTypes.Wye)
%                             %% Wye Connected Voltage Source, Wye Connected Load
%                             for k = 1:nPhases
%                                 I = (ind + 1):(ind + nSeries * nParallel);
%                                 
%                                 rows(I) = reshape(connectionMatrices{j}, 1, []);
%                                 cols(I) = k;
%                                 vals(I) = reshape(connectionPolarity{j}, 1, []) .* (connectionConductivity / assembly.Length.Value / nSeries * assembly.ModeledFraction);
%                                 if j ~= k
%                                     vals(I) = vals(I) * (  - parallelResistance / phaseResistance(k));
%                                 else
%                                     vals(I) = vals(I) * (1 - parallelResistance / phaseResistance(k));
%                                 end
%                                 
%                                 ind = ind + nSeries * nParallel;
%                             end
%                         elseif (sources(i).Type == SourceTypes.VoltageSource) && (assembly.ConnectionType == ConnectionTypes.Delta)
%                             %% Wye Connected Voltage Source, Delta Connected Load
%                           	I = (ind + 1):(ind + nSeries * nParallel);
%                             
%                             rows(I) = reshape(connectionMatrices{j}, 1, []);
%                             cols(I) = j;
%                             vals(I) = reshape(connectionPolarity{j}, 1, []) / phaseResistance(j);
%                             vals(I) = vals(I) .* (connectionConductivity / assembly.Length.Value / nSeries * assembly.ModeledFraction);
%                             
%                             ind = ind + nSeries * nParallel;
%                             
%                             I = (ind + 1):(ind + nSeries * nParallel);
%                             
%                             warning('MotorProto:Verbose', 'This can be generalized to allow all possible m-phase wye-delta connections');
%                             k = mod(j + ceil(nPhases / 2) - 1, nPhases) + 1; %maximum voltage connection
%                             
%                             rows(I) = reshape(connectionMatrices{j}, 1, []);
%                             cols(I) = k;
%                             vals(I) = - reshape(connectionPolarity{j}, 1, []) / phaseResistance(j);
%                             vals(I) = vals(I) .* (connectionConductivity / assembly.Length.Value / nSeries * assembly.ModeledFraction);
%                             
%                             ind = ind + nSeries * nParallel;
%                         else
%                             error('MotorProto:MatrixFactory', 'No implementation for %s connected %s with a %s connected load',...
%                                     char(sources(i).ConnectionType), char(sources(i).Type), char(assembly.ConnectionType));
%                         end
%                     end
%                 elseif (sources(i).ConnectionType == ConnectionTypes.Delta)
%                     error('No Implementation');
%                 else
%                     error('MotorProto:MatrixFactory','No implementation for %s connected sources',char(sources(i).ConnectionType));
%                 end
%                 
%                 rows(ind+1:end) = [];
%                 cols(ind+1:end) = [];
%                 vals(ind+1:end) = [];
%                 
%                 Sr2J = sparse(rows, cols, vals, nComponents, nPhases);
%             end

            if nSources > 0
                mqsInput(iMesh).Sr2El = J2Node * this.PostProcessing(iMesh).F2J;
            else
                mqsInput(iMesh).Sr2El = sparse(nUnknowns,1);
            end
            
            %% UpdateStructure
            structure.Magnetic = mqsInput;
        end
        
        %% Boundary Conditions
        function this = applyBoundaryConditions(this)
            this = applyTangentialBoundaryConditions(this);
            this = applyRadialBoundaryConditions(this);
        end
        
        function mats = applyPeriodicBoundaryConditions(this, mats, iMesh, action)
            if ~strcmpi(action,'none')
                mesh  = this.Mesh_(iMesh);
                I     = mesh.PeriodicBoundaryNodes(1, :);
                J     = mesh.PeriodicBoundaryNodes(2, :);
                N     = numel(mesh.X);
                K     = setdiff(1:N, [I, J]);
                
                if isstruct(mats)
                    field = fields(mats);
                    for i = 1:numel(field)
                        switch lower(action)
                            case {'row', 'rows'}
                                mats(iMesh).(field{i})(I,:) = mats(iMesh).(field{i})(I,:) - mats(iMesh).(field{i})(J,:);
                                mats(iMesh).(field{i})(J,:) = [];
                            case {'column', 'columns'}
                                mats(iMesh).(field{i})(:,I) = mats(iMesh).(field{i})(:,I) - mats(iMesh).(field{i})(:,J);
                                mats(iMesh).(field{i})(:,J) = [];
                            case 'both'
                                mats(iMesh).(field{i})(I,I) = mats(iMesh).(field{i})(I,I) + mats(iMesh).(field{i})(J,J);
                                mats(iMesh).(field{i})(I,K) = mats(iMesh).(field{i})(I,K) - mats(iMesh).(field{i})(J,K);
                                mats(iMesh).(field{i})(K,I) = mats(iMesh).(field{i})(K,I) - mats(iMesh).(field{i})(K,J);
                                mats(iMesh).(field{i})(J,:) = [];
                                mats(iMesh).(field{i})(:,J) = [];
                            case 'reducerow'
                                mats(iMesh).(field{i})(J,:) = [];
                            case 'restorerow'
                                mats(iMesh).(field{i})(J,:) = -mats(iMesh).(field{i})(I,:);
                                mats(iMesh).(field{i})(:,J) = [];
                            case 'boundaryoperation'
                                mats(iMesh).(field{i})      = mats(iMesh).(field{i}) / mesh.Assembly.ModeledFraction;
                                mats(iMesh).(field{i})(:,I) = mats(iMesh).(field{i})(:,I) - mats(iMesh).(field{i})(:,J);
                                mats(iMesh).(field{i})(:,J) = [];
                            otherwise
                                error('MotorProto:MatrixFactory', 'Unknown argument "%s" for method applyPeriodicBoundaryConditions', action);
                        end
                    end
                else
                    switch lower(action)
                        case {'row', 'rows'}
                            mats(I,:) = mats(I,:) - mats(J,:);
                            mats(J,:) = [];
                        case {'column', 'columns'}
                            mats(:,I) = mats(:,I) - mats(:,J);
                            mats(:,J) = [];
                        case 'both'
                            mats(I,I) = mats(I,I) + mats(J,J);
                            mats(I,K) = mats(I,K) - mats(J,K);
                            mats(K,I) = mats(K,I) - mats(K,J);
                            mats(J,:) = [];
                            mats(:,J) = [];
                        case 'reducerow'
                            mats(J,:) = [];
                        case 'restorerow'
                            mats(J,:) = -mats(I,:);
                            mats(:,J) = [];
                        otherwise
                            error('MotorProto:MatrixFactory', 'Unknown argument "%s" for method applyPeriodicBoundaryConditions', action);
                    end
                end
            end
        end
        
        function this = applyTangentialBoundaryConditions(this)
            warning('MotorProto:Verbose', 'When generalizing beyond the assumption of antiperiodic symmetry, check for symmetry type here');
          	%% apply periodic boundary conditions to curl matrices
            %  remember in the three-dimensional case we will need to
            %  perform a rotation of the x and y components of the element
            %  to node curl matrices
            
            mesh  = this.Mesh_;
            index = this.Index;
            nMesh = numel(mesh);
            for i = 1:nMesh
%                 n = numel(mesh(i).X);
%                 I = mesh(i).PeriodicBoundaryNodes(1,:);
                J = mesh(i).PeriodicBoundaryNodes(2,:);
                
            	%% Indices
                warning('MotorProto:Verbose', 'Create general function to remove indices / nodes');
                    
                index.Local(i).A(J)  = [];
                index.Global(i).A(J) = [];
                d = cumsum(diff(index.Local(i).A) - 1);
                d = [0,d] + min(index.Local(i).A) - 1;
                
                index.Local(i).A  = index.Local(i).A  - d;
                index.Global(i).A = index.Global(i).A - d;
                
                if i ~= 1
                    index.Global(i).A = index.Global(i).A - (min(index.Global(i).A) - max(index.Global(i-1).A) - 1);
                end
                
                index.Local(i).X(J)  = [];
                index.Global(i).X(J) = [];
                d = cumsum(diff(index.Local(i).X) - 1);
                d = [0,d] + min(index.Local(i).X) - 1;
                
                index.Local(i).X  = index.Local(i).X  - d;
                index.Global(i).X = index.Global(i).X - d;
                
                if i ~= 1
                    index.Global(i).X = index.Global(i).X - (min(index.Global(i).X) - max(index.Global(i-1).X) - 1);
                end
                
                index.Local(i).Unknowns = index.Local(i).Unknowns  - numel(J);
                
                for j = i:nMesh
                    index.Global(j).Unknowns = index.Global(j).Unknowns - numel(J);
                end
            end
            this.Index = index;
        end
        
        function this = applyRadialBoundaryConditions(this)
            for i = 1:numel(this.Mesh_)
                this.Stiffness.Reluctivity(i).El2El ...
                    = this.Stiffness.Reluctivity(i).El2El ...
                        + real(this.Boundary.Radial.S(i).Inner * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner) ...
                        + real(this.Boundary.Radial.S(i).Outer * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer);
            end
        end
        
      	function harmonics = radialBoundaryHarmonics(this,iMesh)
            warning('Use getRadialBoundaryHarmonics');
            harmonics = getRadialBoundaryHarmonics(this,iMesh);
        end
        
        %% Master Loop Function
      	function structure = buildMatrices(this, structure, methodName)
            for iMesh = numel(this.Mesh_):-1:1;
                [structure, newFields, pbcType] = this.(methodName)(structure, iMesh);
                for iField = 1:numel(newFields)
                    structure.(newFields{iField}) = applyPeriodicBoundaryConditions(this, structure.(newFields{iField}), iMesh, pbcType{iField});
                end
            end
        end
        
        %% Auxillary Functons
        function [t, h]    = getTimePoints(this, Nt)
            warning('MotorProto:Verbose', 'Present implementation assumes synchronous operation');
            [t,h] = this.Assemblies_.getTimePoints(Nt);
        end
        
        function harmonics = getRadialBoundaryHarmonics(this, iMesh)
           	%% Determine the number of boundaries and edges
            edges       = [this.Mesh_.RadialBoundaryEdges];
            nBoundaries = numel(edges);
            nEdges      = cellfun('prodofsize',edges) / 2;
            
            %% Determine the range and spacing of spatial harmonics to include in the Fourier expansion
            nMin = this.Model_.SpatialSymmetries;
            if this.Model_.HasHalfWaveSymmetry
                dn = 2 * nMin;
            else
                dn = nMin;
            end
            nMax = 2 * nEdges * nMin;
            
            %% Calculate harmonics
            harmonics    = cell(1,nBoundaries);
            I            = nMin:dn:nMax(1);
            harmonics{1} = [-fliplr(I),I];
            for i = 2:2:(nBoundaries - 1)
                I              = nMin:dn:max(nMax(i:i+1));
                harmonics{i}   = [-fliplr(I),I];
                harmonics{i+1} = harmonics{i};
            end
            I              = nMin:dn:nMax(end);
            harmonics{end} = [-fliplr(I),I];
            
            if nargin == 2
                harmonics = harmonics((2*iMesh-1):(2*iMesh));
            end
        end
        
     	function gHandleOut = plot(this)
            regions  = this.Regions_;
            nRegions = numel(regions);
            elements = this.Elements;
            triX     = this.X(elements);
            triY     = this.Y(elements);
            elRegion = this.ElementRegions;
            fillArgs = cell(3*nRegions,1);

            for iRegion = 1:nRegions
                elInRegion        	 = (elRegion == iRegion);
                iFillArg             = (iRegion-1)*3;
                fillArgs{iFillArg+1} = triX(:,elInRegion);
                fillArgs{iFillArg+2} = triY(:,elInRegion);
                fillArgs{iFillArg+3} = regions(iRegion).Material.Color;
            end
            
            if nargout == 1
                gHandleOut = fill(fillArgs{:});
            else
                fill(fillArgs{:});
            end
            
            axis equal;
        end
    end
    
    methods (Abstract)
     	%% Matrix and Vector Functions
        linearMatrix       = K(this,t,h)
        
        conductivityMatrix = C(this,t,h)
        
        exogenousFunction  = f(this,t,h)
        
        [nonlinearJacobian, nonlinearFunction] = G(this,t,x)
    end
    
    methods%% Post Processing
       	%% Field Variables
        function [x, x_t, labels, text, nTimes] = continuumVariablePreProcessing(~, solver, dataType, dataPoints)
         	times = solver.Times;
            
            switch lower(dataType)
                case 'time'
                    x      = solver.X(:,dataPoints);
                    x_t    = solver.X_t(:,dataPoints);
                    labels = times(dataPoints);
                    nTimes = numel(dataPoints);
                    text   = 't = %0.4g';
                case 'harmonic'
                	x      = solver.X;
                    x_t    = solver.X_t;
                    labels = dataPoints;
                    nTimes = numel(times);
                    text   = 'h = %d';
                otherwise
                    x      = solver.X;
                    x_t    = solver.X_t;
                    labels = '';
                    nTimes = numel(times);
                    text   = '';
            end
        end
        
        function [a, labels, text] = A(this, solver, dataType, dataPoints)
            [x, ~, labels, text, nTimes] = continuumVariablePreProcessing(this, solver, dataType, dataPoints);
            
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            a           = cell(nAssemblies, nTimes);
            nRows       = zeros(nAssemblies, 1);
            for i = 1:nAssemblies
                I = ppMatrices(i).X2A;
                for j = 1:nTimes
                    a{i,j} = I * x{i,j};
                end
                nRows(i) = length(a{i,1});
            end
            
            if strcmpi(dataType,'harmonic')
                a = cell2mat(a(:,1:(end-1)));
                a = fft(a,[],2) / (nTimes - 1);
                a = mat2cell(a,nRows,ones(1,nTimes-1));
                a = a(:,mod(dataPoints,nTimes - 1) + 1);
            end
        end
        
        function [e, labels, text] = E(this, solver, dataType, dataPoints)
            [x, x_t, labels, text, nTimes] = continuumVariablePreProcessing(this,solver,dataType,dataPoints);
            
            if strcmpi(dataType,'time')
                t = labels;
            else
                t = solver.Times;
            end
            
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            e           = cell(nAssemblies,nTimes);
            assembly    = this.Assemblies_;
            
            nRows       = zeros(nAssemblies,1);
          	nVertPerEl  = numel(ppMatrices(1).X_t2E);
            
            for i = 1:nAssemblies
                for j = 1:nTimes
                    e{i,j} = ppMatrices(i).X2E * x{i,j};

                    sources = assembly(i).Sources;
                    if numel(sources) > 0
                        e{i,j} = e{i,j} + ppMatrices(i).F2E * sources.f(t(j));
                    end

                    e{i,j}     = repmat(e{i,j},1,nVertPerEl);

                    for k = 1:nVertPerEl
                        e{i,j}(:,k) = e{i,j}(:,k) + ppMatrices(i).X_t2E{k} * x_t{i,j};
                    end
                end
                nRows(i) = length(e{i,1});
            end
            
            if strcmpi(dataType,'harmonic')
                e = cellfun(@(x)(reshape(x,[],1)),e,'UniformOutput',false);
                e = cell2mat(e(:,1:(end-1)));
                e = fft(e,[],2) / (nTimes - 1);
                e = mat2cell(e,nRows*nVertPerEl,ones(1,nTimes-1));
                e = e(:,mod(dataPoints,nTimes - 1) + 1);
                e = cellfun(@(x)(reshape(x,[],nVertPerEl)),e,'UniformOutput',false);
            end
        end
        
        function [j, labels, text] = J(this, solver, dataType, dataPoints)
            [x, x_t, labels, text, nTimes] = continuumVariablePreProcessing(this,solver,dataType,dataPoints);
            
            if strcmpi(dataType,'time')
                t = labels;
            else
                t = solver.Times;
            end
            
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            j           = cell(nAssemblies,nTimes);
            assembly    = this.Assemblies_;
            nRows       = zeros(nAssemblies,1);
            nVertPerEl  = numel(ppMatrices(1).X_t2J);
            
            for i = 1:nAssemblies
                for k = 1:nTimes
                    j{i,k} = ppMatrices(i).X2J * x{i,k};

                    sources = assembly(i).Sources;
                    if numel(sources) > 0
                        j{i,k} = j{i,k} + ppMatrices(i).F2J * sources.f(t(k));
                    end

                    j{i,k} = repmat(j{i,k},1,nVertPerEl);

                    for l = 1:nVertPerEl
                        j{i,k}(:,l) = j{i,k}(:,l) + ppMatrices(i).X_t2J{l} * x_t{i,k};
                    end
                end
                nRows(i) = length(j{i,1});
            end

         	if strcmpi(dataType,'harmonic')
                j = cellfun(@(x)(reshape(x,[],1)),j,'UniformOutput',false);
                j = cell2mat(j(:,1:(end-1)));
                j = fft(j,[],2) / (nTimes - 1);
                j = mat2cell(j,nRows*nVertPerEl,ones(1,nTimes-1));
                j = j(:,mod(dataPoints,nTimes - 1) + 1);
                j = cellfun(@(x)(reshape(x,[],nVertPerEl)),j,'UniformOutput',false);
            end
        end
        
        function [b, labels, text] = B(this, solver, dataType, dataPoints)
            [bx, by, labels, text] = BVector(this, solver, dataType, dataPoints);
            
            b = cellfun(@(x,y)(hypot(x,y)),bx,by,'UniformOutput',false);
        end
        
        function [bx, by, labels, text] = BVector(this, solver, dataType, dataPoints)
            [x, ~, labels, text, nTimes] = continuumVariablePreProcessing(this, solver, dataType, dataPoints);
            
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            bx          = cell(nAssemblies,nTimes);
            by          = cell(nAssemblies,nTimes);

            del         = this.Del;
            nRows       = zeros(nAssemblies,1);
            for i = 1:nAssemblies
                for j = 1:nTimes
                    bx{i,j} = del.Curl(i).Zn2Xe * x{i,j};
                    by{i,j} = del.Curl(i).Zn2Ye * x{i,j};
                end
                nRows(i) = length(bx{i,1});
            end
            
         	if strcmpi(dataType,'harmonic')
                bx = cell2mat(bx(:, 1:(end-1)));
                bx = fft(bx, [], 2) / (nTimes - 1);
                bx = mat2cell(bx, nRows, ones(1,nTimes-1));
                bx = bx(:,mod(dataPoints, nTimes - 1) + 1);
                
                by = cell2mat(by(:, 1:(end-1)));
                by = fft(by, [], 2) / (nTimes - 1);
                by = mat2cell(by, nRows, ones(1,nTimes-1));
                by = by(:, mod(dataPoints, nTimes - 1) + 1);
            end
        end
        
       	function [m, labels, text] = M(this, solver, dataType, dataPoints)
            [mx, my, labels, text] = MVector(this, solver, dataType, dataPoints);
            
            m = cellfun(@(x,y)(hypot(x,y)),mx,my,'UniformOutput',false);
        end
        
        function [mx, my, labels, text] = MVector(this, solver, dataType, dataPoints)
            [x, ~, labels, text, nTimes] ...
                = continuumVariablePreProcessing(this,solver,dataType,dataPoints);
            
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            
            mx          = cell(nAssemblies,nTimes);
            my          = cell(nAssemblies,nTimes);
            
            del         = this.Del;
            assembly    = this.Assemblies_;
            mesh        = this.Mesh_;
            nRows       = zeros(nAssemblies,1);
            
            for i = 1:nAssemblies
                for j = 1:nTimes
                    bx = del.Curl(i).Zn2Xe * x{i,j};
                    by = del.Curl(i).Zn2Ye * x{i,j};
                    s  = size(bx);

                    mx{i,j} = zeros(s);
                    my{i,j} = zeros(s);
                    materials = [assembly(i).Regions.Material];
                    for k = 1:numel(materials)
                        K = (mesh(i).ElementRegions == k);
                        [mx{i,j}(K),my{i,j}(K)] = materials(k).vectorM(bx(K),by(K));
                    end
                    [~,rMx,rMy] = materials.remnantMagnetization;
                    mx{i,j}     = mx{i,j} + rMx(mesh(i).ElementRegions);
                    my{i,j}     = my{i,j} + rMy(mesh(i).ElementRegions);
                end
                nRows(i) = length(mx{i,1});
            end
            
         	if strcmpi(dataType,'harmonic')
                mx = cell2mat(mx(:,1:(end-1)));
                mx = fft(mx,[],2) / (nTimes - 1);
                mx = mat2cell(mx,nRows,ones(1,nTimes-1));
                mx = mx(:,mod(dataPoints,nTimes - 1) + 1);
                
                my = cell2mat(my(:,1:(end-1)));
                my = fft(my,[],2) / (nTimes - 1);
                my = mat2cell(my,nRows,ones(1,nTimes-1));
                my = my(:,mod(dataPoints,nTimes - 1) + 1);
            end
        end
        
        function [hx, hy, labels, text] = HVector(this, solver, dataType, dataPoints)
            [bx, by, labels, text] = BVector(this, solver, dataType, dataPoints);
            [mx, my]               = MVector(this, solver, dataType, dataPoints);
            
            hx = cellfun(@(b,m)(b/mu_o - m), bx, mx, 'UniformOutput',false);
            hy = cellfun(@(b,m)(b/mu_o - m), by, my, 'UniformOutput',false);
        end
        
     	function [h, labels, text] = H(this, solver, dataType, dataPoints)
            [hx, hy, labels, text] = HVector(this, solver, dataType, dataPoints);
            
            h = cellfun(@(x, y)(hypot(x, y)), hx, hy, 'UniformOutput', false);
        end
        
        %% Density Variables
        function [l, labels, text] = AverageLossDensity(this, solver, varargin)
            [lcond, lcondtot] = AverageConductionLossDensity(this, solver);
            [lcore, lcoretot] = AverageCoreLossDensity(this, solver);
            
            l = lcond;
            for i = 1:numel(l)
                l{i} = bsxfun(@plus, l{i}, lcore{i});
            end
            
            text   = '(P_{loss} = %0.3g W)';
            labels = lcondtot + lcoretot;
        end
        
        function [l, labels, text] = AverageConductionLossDensity(this, solver, varargin)
            l      = ConductionLossDensity(this, solver, 'Harmonic', 0);
            text   = '(P_{cond} = %0.3g W)';
            labels = AverageConductionLosses(this, solver, 'time');
            labels = labels{1};
        end
        
        function [l, labels, text] = AverageCoreLossDensity(this, solver, varargin)
            pHarmonics = CoreLossDensity(this, solver);
            text       = '(P_{core} = %0.3g W)';
            
            mesh   = this.Mesh_;
            nMesh  = numel(mesh);
            l      = cell(nMesh,1);
            labels = zeros(1,nMesh);
            for i = 1:nMesh
                l{i}      = sum([pHarmonics{i,:}],2);
                labels(i) = mesh(i).ElementAreas * l{i} * this.Model_.Assemblies(i).Length.Value / this.Model_.Assemblies(i).ModeledFraction;
            end
        end
                
      	function [l, labels, text] = ConductionLossDensity(this, solver, dataType, dataPoints)
            warning('MotorProto:Verbose','Need to take into account shape functions when multiplying E and J. Quadratic terms are introduced');
                
            [~, ~, text,labels] = continuumVariablePreProcessing(this, solver, dataType, dataPoints);
            switch lower(dataType)
                case 'time'
                    e = E(this, solver, 'time', dataPoints);
                    j = J(this, solver, 'time', dataPoints);
                case 'harmonic'
                    [~, text, labels] = E(this, solver, 'harmonic', dataPoints);
                    times  = solver.Times;
                    nTimes = numel(times);
                    e = E(this, solver, 'time', 1:nTimes);
                    j = J(this, solver, 'time', 1:nTimes);
            end
            l     = cellfun(@(x,y)(x.*y),e,j,'UniformOutput',false);
            nRows = cellfun('length',l(:,1));
            [~,nVertPerEl] = size(l{1,1});
            
            if strcmpi(dataType,'harmonic')
                l = cellfun(@(x)(reshape(x,[],1)),l,'UniformOutput',false);
                l = cell2mat(l(:,1:(end-1)));
                l = fft(l,[],2) / (nTimes - 1);
                l = mat2cell(l,nRows*nVertPerEl,ones(1,nTimes-1));
                l = l(:,mod(dataPoints,nTimes - 1) + 1);
                l = cellfun(@(x)(reshape(x,[],nVertPerEl)),l,'UniformOutput',false);
            end
        end
        
        function [l, labels, text] = CoreLossDensity(this, solver, varargin)
            times      = solver.Times;
            nHarmonics = numel(times - 1) / 2 - 1;
            h          = 0:(nHarmonics - 1);
            
            [b, text, labels] = this.B(solver, 'Harmonic', h);
            f_e = [this.Assemblies_.ElectricalFrequency];
            f_e = [f_e.Value];
            
            assert( all(f_e - mean(f_e) < sqrt(eps) * mean(f_e)) || all(f_e == 0), 'MotorProto:StaticMatrixFactory', 'All electrical frequencies should be identical');
            
            f_e = mean(f_e);
            f   = f_e * h;
            
            mesh   = this.Mesh_;
            nMesh  = numel(mesh);
            l      = cellfun(@(x)(x*0),b,'UniformOutput',false);
            for i = 1:nMesh
                material   = mesh(i).Materials;
                nMaterials = numel(material);
                for j = 1:nMaterials
                  	J = (mesh(i).ElementRegions == j);
                    if ~material(j).Linear
                        s = material(j).CoreLossCoefficients;
                        for k = 1:nHarmonics
                            l{i,k}(J) = l{i,k}(J)+s(1)*f(k)^s(2)*(2*b{i,k}(J)).^s(3);
                        end
                    else
                        for k = 1:nHarmonics
                            l{i,k}(J) = 0;
                        end
                    end
                end
            end            
        end
        
        %% Bulk Power Variables
        function [p, figLabels, figTitles]   = AverageLosses(this, solver, dataType)
            if nargin < 3
                dataType = [];
            end
            
            p_cond = AverageConductionLosses(this, solver, dataType);
            p_core = AverageCoreLosses(this, solver, dataType);
            p      = sum(p_cond{1}) + sum(p_core{1});
            
            p         = {p};
            figLabels = {''};
            figTitles = {''};
        end

        function [l, figLabels, figTitles]   = AverageCoreLosses(this, solver, dataType)            
            switch lower(dataType)
                case {'default','time'}
                    pHarmonics = CoreLossDensity(this, solver);
                    mesh   = this.Mesh_;
                    nMesh  = numel(mesh);
                    l      = zeros(1,nMesh);
                    for i = 1:nMesh
                        l(i) = mesh(i).ElementAreas * sum([pHarmonics{i,:}],2) * this.Model_.Assemblies(i).Length.Value / this.Model_.Assemblies(i).ModeledFraction;
                    end
            
                    l         = {l};
                    figLabels = {''};
                    figTitles = {''};
                case 'harmonic'
                    error('No Implementation');
            end
        end
        
     	function [l, figLabels, figTitles]   = InstantaneousConductionLosses(this, solver, dataType)
            x           = solver.X;
            x_t         = solver.X_t;
            t           = solver.Times;
            Nt          = numel(t);
            mesh        = this.Mesh_;
            assembly    = this.Assemblies_;
            nAssemblies = numel(assembly);
            l           = zeros(nAssemblies,numel(t));
            ppMatrix    = this.PostProcessing;
            
            for p = 1:nAssemblies;
                source   = [assembly(p).Sources];
                nSources = numel(source);
                ela      = mesh(p).ElementAreas.';
                for k = 1:Nt;
                    e0 = ppMatrix(p).X2E * x{p,k};
                    j0 = ppMatrix(p).X2J * x{p,k};
                    
                    if nSources > 0
                        e0 = e0 + ppMatrix(p).F2E * source.f(t(k));
                        j0 = j0 + ppMatrix(p).F2J * source.f(t(k));
                    end
                    
                    ne = numel(ppMatrix(p).X_t2E);
                    nj = numel(ppMatrix(p).X_t2J);

                    e1 = zeros(size(e0));
                    for i = 1:ne
                        e1(:,i) = ppMatrix(p).X_t2E{i} * x_t{p,k};
                    end

                    j1 = zeros(size(j0));
                    for i = 1:nj
                        j1(:,i) = ppMatrix(p).X_t2J{i} * x_t{p,k};
                    end

                    %% 0 x 0
                    l(p,k) = j0'*(ela.*e0);

                    %% 1 x 0, 0 x 1
                    for i = 1:ne
                        c = 2 * factorial(0) * factorial(1) / factorial(2 + 0 + 1);
                        l(p,k) = l(p,k) + c*j0'*(ela.*e1(:,i));
                    end

                    for i = 1:nj
                        c = 2 * factorial(0) * factorial(1) / factorial(2 + 0 + 1);
                        l(p,k) = l(p,k) + c*j1(:,i)'*(ela.*e0);
                    end

                    %% 1 x 1
                    for i = 1:ne
                        for j = 1:ne
                            if i == j
                                ik = 2;
                                jk = 0;
                            else
                                ik = 1;
                                jk = 1;
                            end
                            c = 2*factorial(ik)*factorial(jk)/factorial(2+ik+jk);
                            l(p,k) = l(p,k) + c*j1(:,j)'*(ela.*e1(:,i));
                        end
                    end
                end
                
                l(p,:) = l(p,:) * this.Model_.Assemblies(p).Length.Value / this.Model_.Assemblies(p).ModeledFraction;
            end
            
            if strcmpi(dataType,'harmonic')
                l = fft(l(:,1:(Nt-1)),[],2) / (Nt - 1);
            end
            
            l         = num2cell(l,2);
            figLabels = {''};
            figTitles = {''};
        end
        
        function [l, figLabels, figTitles]   = AverageConductionLosses(this, solver, dataType)
            switch lower(dataType)
                case {'default','time'}
                    [l, figLabels, figTitles] = InstantaneousConductionLosses(this, solver, 'Time');
                    l = cellfun(@(x)(mean(x)), l).';
                    l = {l};
                case 'harmonic'
                   	x        = solver.X;
                    x_t      = solver.X_t;
                    t        = solver.Times;
                    Nt       = numel(t);
                    nRows    = cellfun('length', x);
                    nRows    = nRows(:,1);
                    
                    x        = cell2mat(x);
                    x        = fft(full(x(:,1:(Nt-1))),[],2) / (Nt - 1);
                    x        = mat2cell(x,nRows,ones(1,Nt-1));
                    
                    x_t      = cell2mat(x_t);
                    x_t      = fft(full(x_t(:,1:(Nt-1))),[],2) / (Nt - 1);
                    x_t      = mat2cell(x_t,nRows,ones(1,Nt-1));
                    
                    mesh     = this.Mesh_;
                    assembly = this.Assemblies_;
                    l        = zeros(1,numel(t));
                    ppMatrix = this.PostProcessing;

                    for p = 1:numel(ppMatrix);
                        source   = [assembly(p).Sources];
                        f        = source.f(t);
                        f        = fft(f(:,1:(Nt-1)),[],2) / (Nt - 1);
                        nSources = numel(source);
                        ela      = mesh(p).ElementAreas.';
                        for k = 1:(Nt-1);
                            e0 = ppMatrix(p).X2E * x{p,k};
                            j0 = ppMatrix(p).X2J * x{p,k};

                            if nSources > 0
                                e0 = e0 + ppMatrix(p).F2E * f(:,k);
                                j0 = j0 + ppMatrix(p).F2J * f(:,k);
                            end

                            ne = numel(ppMatrix(p).X_t2E);
                            nj = numel(ppMatrix(p).X_t2J);

                            e1 = zeros(size(e0));
                            for i = 1:ne
                                e1(:,i) = ppMatrix(p).X_t2E{i} * x_t{p,k};
                            end

                            j1 = zeros(size(j0));
                            for i = 1:nj
                                j1(:,i) = ppMatrix(p).X_t2J{i} * x_t{p,k};
                            end

                            %% 0 x 0
                            l(k) = l(k) + j0'*(ela.*e0);

                            %% 1 x 0, 0 x 1
                            for i = 1:ne
                                c = 2 * factorial(0) * factorial(1) / factorial(2 + 0 + 1);
                                l(k) = l(k) + c*j0'*(ela.*e1(:,i));
                            end

                            for i = 1:nj
                                c = 2 * factorial(0) * factorial(1) / factorial(2 + 0 + 1);
                                l(k) = l(k) + c*j1(:,i)'*(ela.*e0);
                            end

                            %% 1 x 1
                            for i = 1:ne
                                for j = 1:ne
                                    if i == j
                                        ik = 2;
                                        jk = 0;
                                    else
                                        ik = 1;
                                        jk = 1;
                                    end
                                    c = 2*factorial(ik)*factorial(jk)/factorial(2+ik+jk);
                                    l(k) = l(k) + c*j1(:,j)'*(ela.*e1(:,i));
                                end
                            end
                        	l(k) = l(k) * this.Model_.Assemblies(1).Length.Value / this.Model_.SpaceModelFraction;
                        end
                    end
                 	l         = {l};
                    figLabels = {''};
                    figTitles = {''};
            end    
        end
        
        function [tau, figLabels, figTitles] = Torque(this, solver, dataType)
        	warning('MotorProto:Verbose', 'Generalize length calculation');
            warning('MotorProto:Verbose', 'Generalize torque calculation to multiple annulii');
            
            x    = solver.X;
            t    = solver.Times;
            ppm  = this.PostProcessing;
            
            R    = this.Boundary.Radial.R;
            D    = this.Boundary.Radial.D;
            F    = this.Boundary.Radial.F;
            G    = this.Boundary.Radial.G;
            
            tau  = zeros(2,length(t));
            nh   = this.getRadialBoundaryHarmonics(2);
            nh   = nh{1}.';
            r    = this.Assemblies_(2).InnerRadius.Value;
            l    = this.Assemblies_(2).Length.Value;
            Nt   = length(t);
                
            if iscell(D(1).Outer)
                [xRe, xIm] = dscfft(cell2mat(x(1,1:(end-1))),[],2);
                xRe = ppm(1).Full2Reduced * xRe;
                xIm = ppm(1).Full2Reduced * xIm;
                hRe = cell(1,(Nt - 1) / 2);
                hIm = cell(1,(Nt - 1) / 2);
                for i = 1:((Nt - 1) / 2)
                    j      = mod(i-1,numel(D(1).Outer)) + 1;
                    hRe{i} = D(1).Outer{j}{1,1} * xRe(:,i) + D(1).Outer{j}{1,2} * xIm(:,i);
                    hIm{i} = D(1).Outer{j}{2,1} * xRe(:,i) + D(1).Outer{j}{2,2} * xIm(:,i);
                end
                hRe     = cell2mat(hRe);
                hIm     = cell2mat(hIm);
                htOuter = [hRe + 1i*hIm, hIm(:,1), fliplr(hRe(:,2:end) - 1i*hIm(:,2:end))];
                htOuter = ifft(htOuter, [], 2) * (Nt-1);
                
                [xRe, xIm] = dscfft(cell2mat(x(2,1:(end-1))),[],2);
                xRe = ppm(2).Full2Reduced * xRe;
                xIm = ppm(2).Full2Reduced * xIm;
                hRe = cell(1,(Nt - 1) / 2);
                hIm = cell(1,(Nt - 1) / 2);
                for i = 1:((Nt - 1) / 2)
                    j      = mod(i-1,numel(D(2).Inner)) + 1;
                    hRe{i} = D(2).Inner{j}{1,1} * xRe(:,i) + D(2).Inner{j}{1,2} * xIm(:,i);
                    hIm{i} = D(2).Inner{j}{2,1} * xRe(:,i) + D(2).Inner{j}{2,2} * xIm(:,i);
                end
                hRe     = cell2mat(hRe);
                hIm     = cell2mat(hIm);
                htInner = [hRe + 1i*hIm, hIm(:,1), fliplr(hRe(:,2:end) - 1i*hIm(:,2:end))];
                htInner = ifft(htInner, [], 2) * (Nt-1);
                
                for i = 1:(Nt-1)
                    ht = exp(1i * (R(2).Inner - R(1).Outer) * t(i)) .* (G(1).Outer * htOuter(:,i)) + (F(2).Inner * htInner(:,i));

                    br = 1i * nh .* htInner(:,i) / r;

                    tau(1,i) = 2*pi*l*r^2*real(sum(conj(ht).*br));
                    tau(2,i) = tau(1,i);
                end
                tau(:,end) = tau(:,1);
            else
                for i = 1:Nt;
                    ht =   exp(1i * (R(2).Inner - R(1).Outer) * t(i)).* ...
                                (G(1).Outer * D(1).Outer * ppm(1).Full2Reduced * x{1,i})...
                         + (F(2).Inner * D(2).Inner * ppm(2).Full2Reduced * x{2,i});

                    br = 1i * nh .* (D(2).Inner * ppm(2).Full2Reduced *  x{2,i}) / r;

                    tau(1,i) = 2*pi*l*r^2*real(sum(conj(ht).*br));
                    tau(2,i) = tau(1,i);
                end
            end
            tau = mean(tau);
            if strcmpi(dataType,'harmonic')
                tau = fft(tau(1:(Nt-1)),[],2) / (Nt - 1);
            end
            tau       = {tau};
            figLabels = {''};
            figTitles = {''};
        end
        
        function [tau, figLabels, figTitles] = AverageTorque(this, solver, dataType)
        	warning('MotorProto:Verbose', 'Generalize this calculation');
            switch lower(dataType)
                case {'default','time'}
                    [tau, figLabels, figTitles] = Torque(this, solver, 'Time');
                    tau                         = {cellfun(@(x)(mean(x(1:(end-1)))), tau)};
                case 'harmonic'
                    error('No Implementation');
            end
        end
        
        function [pow, figLabels, figTitles] = AverageOutputPower(this, solver, dataType)
        	warning('MotorProto:Verbose', 'Generalize this calculation');
            switch lower(dataType)
                case {'default','time'}
                    [tau, figLabels, figTitles] = AverageTorque(this, solver, 'Time');
                    pow                         = {tau{1} * this.Assemblies_(1).AngularVelocity};
                case 'harmonic'
            end
        end
        
        function [eff, figLabels, figTitles] = Efficiency(this, solver, dataType)
        	warning('MotorProto:Verbose', 'Generalize this calculation');
            if nargin < 3
                dataType = [];
            end
            
            p_loss = AverageLosses(this, solver, dataType);
            p_out  = AverageOutputPower(this, solver, dataType);
            
            eff       = {p_out{1} / (p_out{1} + p_loss{1})};
            figLabels = {''};
            figTitles = {''};
        end
        
        %% Other Bulk Variables
        function [v, figLabels, figTitles]      = Voltage(this, solver, dataType)
            x           = solver.X;
            x_t         = solver.X_t;
            t           = solver.Times;
            Nt          = numel(t);
            assembly    = this.Assemblies_;
            nAssemblies = numel(assembly);
            v           = cell(nAssemblies, 1);
            figLabels   = cell(nAssemblies, 1);
            figTitles   = cell(nAssemblies, 1);
            remove      = false(nAssemblies, 1);
            ppMatrices  = this.PostProcessing;
            
            for i = 1:nAssemblies                
                source    = assembly(i).Sources;
                hasSource = (numel(source) > 0);
                if hasSource
                    figTitles{i} = [source.Name ' Phase'];
                    
                    nPhases      = source.Phases.Value;
                    figLabels{i} = cell(1,nPhases);
                    for j = 1:nPhases
                        figLabels{i}{j} = ['Phase ' char(num2str('A')+j-1)];
                    end
                    
                    v{i} = ppMatrices(i).X2V * x{i,1};
                    v{i} = v{i} + ppMatrices(i).X_t2V * x_t{i,1};
                	v{i} = v{i} + ppMatrices(i).F2V * source.f(t(1));
                    v{i} = repmat(v{i},1,Nt);
                    for j = 2:Nt
                        v{i}(:,j) = ppMatrices(i).X2V * x{i,j};
                        v{i}(:,j) = v{i}(:,j) + ppMatrices(i).X_t2V * x_t{i,j}; 
                      	v{i}(:,j) = v{i}(:,j) + ppMatrices(i).F2V * source.f(t(j));
                    end
                
                    if strcmpi(dataType,'harmonic')
                        v{i} = fft(v{i}(:,1:(end-1)), [] ,2) / (Nt - 1);
                    end
                    v{i}         = num2cell(v{i},2);
                else
                    remove(i) = true;
                end
            end
            v(remove)         = [];
            figLabels(remove) = [];
            figTitles(remove) = [];
        end
        
        function [lambda, figLabels, figTitles] = FluxLinkage(this, solver, dataType)
            x           = solver.X;
            t           = solver.Times;
            Nt          = numel(t);
            assembly    = this.Assemblies_;
            nAssemblies = numel(assembly);
            lambda      = cell(nAssemblies, 1);
            figLabels   = cell(nAssemblies, 1);
            figTitles   = cell(nAssemblies, 1);
            remove      = false(nAssemblies, 1);
            ppMatrices  = this.PostProcessing;
            
            for i = 1:nAssemblies                
                source    = assembly(i).Sources;
                hasSource = (numel(source) > 0);
                if hasSource
                    figTitles{i} = [source.Name];
                    
                    nPhases      = source.Phases.Value;
                    figLabels{i} = cell(1,nPhases);
                    for j = 1:nPhases
                        figLabels{i}{j} = ['Phase ' char(num2str('A')+j-1)];
                    end
                    
                    lambda{i} = ppMatrices(i).X2Lambda * x{i,1};
                    lambda{i} = repmat(lambda{i},1,Nt);
                    for j = 2:Nt
                        lambda{i}(:,j) = ppMatrices(i).X2Lambda * x{i,j};
                    end
                
                    if strcmpi(dataType,'harmonic')
                        lambda{i} = fft(lambda{i}(:,1:(end-1)), [] ,2) / (Nt - 1);
                    end
                    lambda{i} = num2cell(lambda{i},2);
                else
                    remove(i) = true;
                end
            end
            lambda(remove)    = [];
            figLabels(remove) = [];
            figTitles(remove) = [];
        end
        
        function [i, figLabels, figTitles]      = Current(this, solver, dataType)
            x           = solver.X;
            x_t         = solver.X_t;
            t           = solver.Times;
            Nt          = numel(t);
            assembly    = this.Assemblies_;
            nAssemblies = numel(assembly);
            i           = cell(nAssemblies, 1);
            figLabels   = cell(nAssemblies, 1);
            figTitles   = cell(nAssemblies, 1);
            remove      = false(nAssemblies, 1);
            ppMatrices  = this.PostProcessing;
            
            for j = 1:nAssemblies
                source    = assembly(j).Sources;
                hasSource = (numel(source) > 0);
                if hasSource
                    figTitles{j} = [source.Name ' Phase'];
                    
                    nPhases      = source.Phases.Value;
                    figLabels{j} = cell(1,nPhases);
                    for k = 1:nPhases
                        figLabels{j}{k} = ['Phase ' char(num2str('A')+k-1)];
                    end
                    
                    i{j} = ppMatrices(j).X2I * x{j,1};
                    i{j} = i{j} + ppMatrices(j).X_t2I * x_t{j,1};
                    i{j} = i{j} + ppMatrices(j).F2I * source.f(t(1));
                    i{j} = repmat(i{j},1,Nt);
                    for k = 2:Nt
                        i{j}(:,k) = ppMatrices(j).X2I * x{j,k};
                        i{j}(:,k) = i{j}(:,k) + ppMatrices(j).X_t2I * x_t{j,k};
                        i{j}(:,k) = i{j}(:,k) + ppMatrices(j).F2I * source.f(t(k));
                    end
                    
                    if strcmpi(dataType,'harmonic')
                        i{j} = fft(i{j}(:,1:(end-1)), [] ,2) / (Nt - 1);
                    end
                    i{j} = num2cell(i{j},2);
                else
                    remove(j) = true;
                end
            end
            i(remove)         = [];
            figLabels(remove) = [];
            figTitles(remove) = [];
        end
    end
    
 	methods (Sealed)
        function copyOut = copy(this)
            nThis   = numel(this);
            copyOut = this;
            for i = 1:nThis
                copyOut(i).Model_ = copy(this.Model_);
            end
        end
    end
end