classdef DynamicMatrixFactory < StaticMatrixFactory
	methods
        %% Constructor
        function this = DynamicMatrixFactory(varargin)
            warning('MotorProto:Verbose','Present implementation assumes anti-periodic symmetry');
            if nargin > 0
                this.Model = varargin{1};
                this = build(this,'space');
            end
        end

        %% Preprocessing Functions
        function [structure, newFields, pbcType] = buildLocalIndexVectors(this, structure, iMesh)
            [structure, newFields, pbcType] = buildLocalIndexVectors@StaticMatrixFactory(this, structure, iMesh);            
            
            local     = structure.Local;
            assembly  = this.Assemblies_(iMesh);
            nUnknowns = local(iMesh).Unknowns;
            
            %% Regions
            regions   = assembly.Regions;
            nRegions  = numel(regions);
            dynamics  = [regions.Dynamics];
            isDynamic = (dynamics == DynamicsTypes.Dynamic);
            nDynamic  = sum(isDynamic);
            
            local(iMesh).Regions            = zeros(1, nRegions);
            local(iMesh).Regions(isDynamic) = (1:nDynamic) + nUnknowns;
                
            %% Components
            components       = assembly.Components;
            nComponents      = numel(components);
            region2Component = assembly.convertIndex('Regions', 1:nRegions, 'Components');
            
            local(iMesh).Components                              = zeros(1, nComponents);
            local(iMesh).Components(region2Component(isDynamic)) = (1:nDynamic) + nUnknowns;
            
            nUnknowns = nUnknowns + nDynamic;
            
            %% Circuits
            sources    = assembly.Sources;
            nSources   = numel(sources);
            if nSources > 0
                [rows,~,~] = sources.circuitMatrix(assembly.ConnectionType);

                local(iMesh).Circuits = (1:max(rows)) + nUnknowns;

                nUnknowns = nUnknowns + max(rows);
            else
                local(iMesh).Circuits = [];
            end
            
            %% Unknowns
            local(iMesh).Unknowns = nUnknowns;
            
            structure.Local = local;
        end        
        
        function [structure, newFields, pbcType] = buildReluctivityMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildReluctivityMatrices@StaticMatrixFactory(this, structure, iMesh);
            
            reluctivity = structure.Reluctivity;
            
            mesh        = this.Mesh_(iMesh);
            assembly    = this.Assemblies_(iMesh);
            regions     = assembly.Regions;
            nRegions    = numel(regions);
            dynamics    = [regions.Dynamics];
            isDynamic   = (dynamics == DynamicsTypes.Dynamic);
            nDynamic    = sum(isDynamic);
            sources     = assembly.Sources;
            nSources    = numel(sources);
            materials   = [regions.Material];
            ind         = this.Index.Local(iMesh);
            
            el          = mesh.Elements;
            elArea      = mesh.ElementAreas;
            elRegions   = mesh.ElementRegions;
            
            elCond    = [materials.sigma];
            elCond    = elCond(elRegions);
            
            nUnknowns = ind.Unknowns;
            
            %% Region to element coupling (through scalar potential)
            i = [el(1,:), el(2,:), el(3,:)];
            j = repmat(ind.Regions(elRegions), 1 , 3);
            s = repmat(elCond.*elArea, 1, 3) / 3;
            
            isZero    = (j == 0);
            i(isZero) = [];
            j(isZero) = [];
            s(isZero) = [];
            
            reluctivity(iMesh).Re2El = sparse(i, j, s, nUnknowns, nUnknowns);
            
            %% Region to region coupling
            %	Use DISCRETE region areas
            s        = zeros(1, nDynamic);
            j        = 1;
            for i = 1:nRegions
                if isDynamic(i)
                    I    = (elRegions == i);
                    s(j) = sum(elArea(I));
                    j    = j + 1;
                end
            end
            
            i = ind.Regions(isDynamic);
            s = s.*[materials(isDynamic).sigma];
            
            % Replace zero values with mean values to avoid poorly conditioned and singular matrices
            isZero    = (s == 0);
            s(isZero) = mean(s(~isZero));
            
            reluctivity(iMesh).Re2Re = sparse(i, i, s, nUnknowns, nUnknowns);
            
            %% External Circuit and Source Coupling
            if nSources > 0
                [rows, cols, vals]       = sources.circuitMatrix(assembly.ConnectionType);
                rows                     = ind.Circuits(rows);
                cols                     = ind.Circuits(cols);
                vals                     = vals / assembly.Length.Value * assembly.ModeledFraction;
                reluctivity(iMesh).Cr2Cr = sparse(rows, cols, vals, nUnknowns, nUnknowns);
                
                [rows, cols, vals]       = sources.scalarPotential2CircuitMatrix;
                rows                     = ind.Circuits(rows);
                cols                     = ind.Components(cols);
                
                reluctivity(iMesh).Re2Cr = sparse(rows, cols, vals, nUnknowns, nUnknowns);
                reluctivity(iMesh).Cr2Re = sparse(cols, rows, vals, nUnknowns, nUnknowns);
            else
                reluctivity(iMesh).Cr2Cr = sparse(nUnknowns, nUnknowns);
                reluctivity(iMesh).Re2Cr = sparse(nUnknowns, nUnknowns);
                reluctivity(iMesh).Cr2Re = sparse(nUnknowns, nUnknowns);
            end
            
            %% make empty matrices
            reluctivity(iMesh).El2Re = sparse(nUnknowns, nUnknowns);
            
            structure.Reluctivity = reluctivity;
        end
        
        function [structure, newFields, pbcType] = buildConductivityMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildConductivityMatrices@StaticMatrixFactory(this, structure, iMesh);
            
            conductivity = structure.Conductivity;
            
            %% Setup
            mesh      = this.Mesh_(iMesh);
            assembly  = this.Assemblies_(iMesh);
            regions   = assembly.Regions;
            dynamics  = [regions.Dynamics];
            isDynamic = (dynamics == DynamicsTypes.Dynamic);
            materials = [regions.Material];
            ind       = this.Index.Local(iMesh);
            
            el        = mesh.Elements;
            elArea    = mesh.ElementAreas;
            elRegions = mesh.ElementRegions;
            elCond    = [materials.sigma];
            elCond    = elCond(elRegions);
            nUnknowns = ind.Unknowns;
            
            %% Ignore elements which have static dynamics
            I         = ismember(elRegions, find(isDynamic));
            el        = el(:, I);
            elArea    = elArea(I);
            elRegions = elRegions(I);
            elRegions = ind.Regions(elRegions);
            elCond    = elCond(I);            

            %% Element to region coupling (spatial integral of A_t)
            i = repmat(elRegions, 1 , 3);
            j = [el(1,:), el(2,:), el(3,:)];
            s = repmat(elCond.*elArea, 1, 3) / 3;
            
            conductivity(iMesh).El2Re = sparse(i, j, s, nUnknowns, nUnknowns);
            
            %% make empty matrices
            conductivity(iMesh).Re2El = sparse(nUnknowns, nUnknowns);
            
            warning('MotorProto:Verbose','Generalize Region to Region Coupling');
            conductivity(iMesh).Re2Re = sparse(nUnknowns, nUnknowns);
            
            structure.Conductivity = conductivity;
        end
        
     	function [structure, newFields, pbcType] = buildMagneticInputMatrices(this, structure, iMesh)
            newFields = {'Magnetic'};
            pbcType   = {'Row'};
            
            if isfield(structure, newFields{1})
                mqsInput = structure.Magnetic;
            else
                mqsInput = struct([]);
            end
            
            mqsInput(iMesh).Sr2El = buildElementInputMatrices(this, iMesh);
            mqsInput(iMesh).Sr2Re = buildRegionInputMatrices(this , iMesh);
            mqsInput(iMesh).Sr2Cr = buildSourceInputMatrices(this , iMesh);
            
            structure.Magnetic = mqsInput;
        end
        
        function sr2el = buildElementInputMatrices(this, iMesh)
            nUnknowns = this.Index.Local(iMesh).Unknowns;
            sources   = this.Assemblies_(iMesh).Sources;
            nSources  = numel(sources);
            
            if nSources > 0
                sr2el = sparse(nUnknowns, sources.Phases.Value);
            else
                sr2el = sparse(nUnknowns, 0);
            end
        end
        
        function sr2re = buildRegionInputMatrices(this, iMesh)
            nUnknowns = this.Index.Local(iMesh).Unknowns;
            sources   = this.Assemblies_(iMesh).Sources;
            nSources  = numel(sources);
            
            if nSources > 0
                sr2re = sparse(nUnknowns, sources.Phases.Value);
            else
                sr2re = sparse(nUnknowns, 0);
            end
        end
        
        function sr2cr = buildSourceInputMatrices(this, iMesh)
            ind       = this.Index.Local(iMesh);
            nUnknowns = ind.Unknowns;
            assembly  = this.Assemblies_(iMesh);
            sources   = assembly.Sources;
            nSources  = numel(sources);
            
            if nSources > 0
                [rows, cols, vals] = sources.source2CircuitMatrix(assembly.ConnectionType);
                rows  = ind.Circuits(rows);
                vals  = vals / assembly.Length.Value * assembly.ModeledFraction;
                sr2cr = sparse(rows, cols, vals, nUnknowns, sources.Phases.Value);
            else
                sr2cr = sparse(nUnknowns, 0);
            end
        end
        
       	function this  = buildPostProcessingMatrices(this)
            warning('MotorProto:Verbose', 'Assumes anti-periodic symmetry');
            
            mesh           = this.Mesh_;
            nMesh          = numel(mesh);
            index          = this.Index;
            postProcessing = struct.empty(0,nMesh);
            
            for i = 1:nMesh
                nUnknowns = index.Local(i).Unknowns;
                
                assembly   = mesh(i).Assembly;
                components = assembly.Components; 
                regions    = assembly.Regions;
                materials  = [regions.Material];        
                elements   = mesh(i).Elements;
                elRegion   = mesh(i).ElementRegions;
                elAreas    = mesh(i).ElementAreas;
                
                        
                sources   = assembly.Sources;
                
                nSources  = length(sources);   
                nEls      = length(elements);
                
                I         = speye(nUnknowns,nUnknowns);
            
              	postProcessing(i).Full2Reduced = this.applyPeriodicBoundaryConditions(I, i, 'ReduceRow');
                postProcessing(i).Reduced2Full = this.applyPeriodicBoundaryConditions(I, i, 'RestoreRow');
                
                %% Solution to Magnetic Vector Potential
                postProcessing(i).X2A = I(index.Local(i).A,:);
                
                %% Electrical Field / Current Density calculations
                rows = 1:nEls;
                cols = elements;
                
                jVal = [materials.sigma];
                jVal = - jVal(elRegion);
                
                eVal            = -ones(1,nEls);
                eVal(jVal == 0) = 0;
                
                x_t2e = cell(1,3);
                x_t2j = cell(1,3);
                for j = 1:3
                    x_t2e{j} = sparse(rows,cols(j,:), eVal, nEls, nUnknowns);
                    x_t2j{j} = sparse(rows,cols(j,:), jVal, nEls, nUnknowns);
                end
                
                postProcessing(i).X_t2E = x_t2e;
                postProcessing(i).X_t2J = x_t2j;
                
                cols = index.Local(i).Regions(elRegion);
                rows = rows(cols ~= 0);
                eVal = eVal(cols ~= 0);
                jVal = jVal(cols ~= 0);
                cols = cols(cols ~= 0);

                x2e = sparse(rows, cols, eVal, nEls, nUnknowns);
                x2j = sparse(rows, cols, jVal, nEls, nUnknowns);

                postProcessing(i).X2E = x2e;
                postProcessing(i).X2J = x2j;

                if nSources > 0
                    postProcessing(i).F2E = sparse(nEls, sources.Phases.Value);
                    postProcessing(i).F2J = sparse(nEls, sources.Phases.Value);
                    postProcessing(i).F2V = sparse(sources.Phases.Value, sources.Phases.Value);
                    postProcessing(i).F2I = sparse(sources.Phases.Value, sources.Phases.Value);
                else
                    postProcessing(i).F2E = sparse(nEls, 0);
                    postProcessing(i).F2J = sparse(nEls, 0);
                    postProcessing(i).F2V = sparse(0, 0);
                    postProcessing(i).F2I = sparse(0, 0);
                end

                if nSources > 0
                    nPhases = sources.Phases.Value;
                    %   Voltage, Current from grad(Phi)
                    rowsP  = zeros(1,0);
                    colsP  = zeros(1,0);
                    valsPV = zeros(1,0);
                    valsPI = zeros(1,0);
                    
                    %	Flux Linkage, Current from A, dAdT
                    rowsA  = zeros(1,0);
                    colsA  = zeros(1,0);
                    valsAF = zeros(1,0); %Flux Linkage
                    valsAI = zeros(1,0); %Current
                    
                    connectionMatrices = sources.ConnectionMatrices;
                    connectionPolarity = sources.ConnectionPolarity;
                    indA = [0 0];
                    indP = [0 0];
                    for j = 1:nPhases
                        connectionComponents = components(connectionMatrices{j});
                        [nSeries, nParallel] = size(connectionMatrices{j});
                        
                        if ~all(size(connectionComponents) == size(connectionMatrices{j}))
                            connectionComponents = connectionComponents.';
                        end
                        
                        for k = 1:nSeries
                            indA(1) = indA(2);
                            indP(1) = indP(2);
                            totalArea = 0;
                            for l = 1:nParallel
                                m         = find(connectionComponents(k,l) == regions);
                                K         = (elRegion == m);
                                N         = sum(K);
                                areas     = repmat(elAreas(K),3,1);
                                areas     = reshape(areas,1,[]);
                                nodes     = elements(:,K);
                                nodes     = reshape(nodes,1,[]);
                                totalArea = totalArea + sum(elAreas(K));
                                
                                %% grad(Phi)
                                rowsP(indP(2)+1)  = j;
                                colsP(indP(2)+1)  = index.Local(i).Components(connectionMatrices{j}(k,l));
                                valsPV(indP(2)+1) = connectionPolarity{j}(k, l) * assembly.Length.Value / assembly.ModeledFraction;
                                valsPI(indP(2)+1) = connectionPolarity{j}(k, l) * sum(elAreas(K)) * materials(m).sigma;
                                indP(2)           = indP(2) + 1;
                                
                                %% A, dA/dt
                                I         = indA(2)+1:indA(2)+3*N;
                                rowsA(I)  = j;
                                colsA(I)  = nodes;
                                valsAF(I) = connectionPolarity{j}(k,l) * areas * assembly.Length.Value / (3 * assembly.ModeledFraction);
                                valsAI(I) = connectionPolarity{j}(k,l) * areas * materials(m).sigma / 3;
                                indA(2)   = indA(2) + 3*N;
                            end
                            valsAF(indA(1)+1:indA(2)) = valsAF(indA(1)+1:indA(2)) / totalArea;
                            valsAI(indA(1)+1:indA(2)) = valsAI(indA(1)+1:indA(2)) / nSeries;
                            valsPV(indP(1)+1:indP(2)) = valsPV(indP(1)+1:indP(2)) / nParallel;
                            valsPI(indP(1)+1:indP(2)) = valsPI(indP(1)+1:indP(2)) / nSeries;
                        end
                    end
                    postProcessing(i).X2V      = sparse(rowsP, colsP, valsPV, nPhases, nUnknowns);
                    postProcessing(i).X2I      = sparse(rowsP, colsP, valsPI, nPhases, nUnknowns);
                    postProcessing(i).X2Lambda = sparse(rowsA, colsA, valsAF, nPhases, nUnknowns);
                    postProcessing(i).X_t2V    = sparse(nPhases, nUnknowns);
                    postProcessing(i).X_t2I    = sparse(rowsA, colsA, valsAI, nPhases, nUnknowns);
                else
                    %   Voltage
                    postProcessing(i).X2V   = sparse(0, nUnknowns);
                    postProcessing(i).X_t2V = sparse(0, nUnknowns);
                    
                    %   Current
                    postProcessing(i).X2I   = sparse(0, nUnknowns);
                    postProcessing(i).X_t2I = sparse(0, nUnknowns);
                    
                    %   Flux Linkage
                    postProcessing(i).X2Lambda = sparse(0, nUnknowns);
                end
            end
            this.PostProcessing = postProcessing;
        end
        
        %% Boundary Conditions
        function this  = applyTangentialBoundaryConditions(this)
            this = applyTangentialBoundaryConditions@StaticMatrixFactory(this);            
            
           	mesh  = this.Mesh_;
            index = this.Index;
            nMesh = numel(mesh);
            for i = 1:nMesh
                J = mesh(i).PeriodicBoundaryNodes(2,:);
                
                %% Indices
                index.Local(i).Regions  = index.Local(i).Regions  - numel(J);
                index.Global(i).Regions = index.Global(i).Regions - numel(J);
            end
            this.Index = index;
        end
        
        %% Matrix/Vector Functions
        function linearMatrix       = K(this, t, h)
            %% create stiffness matrix
            %  periodic boundary conditions are included in the curl matrices
            assembly     = this.Assemblies_;
            nAssemblies  = numel(assembly);
            linearMatrix = cell(nAssemblies,nAssemblies);
            
            for i = 1:nAssemblies
                linearMatrix{i,i} =   this.Stiffness.Reluctivity(i).El2El ...
                                    + this.Stiffness.Reluctivity(i).Re2El ...
                                    + this.Stiffness.Reluctivity(i).El2Re * h ...
                                    + this.Stiffness.Reluctivity(i).Re2Re * h ...
                                    + this.Stiffness.Reluctivity(i).Re2Cr * h ...
                                    + this.Stiffness.Reluctivity(i).Cr2Re * h ...
                                    + this.Stiffness.Reluctivity(i).Cr2Cr * h;
                
                if i > 1
                    R = this.Boundary.Radial.R(i).Inner - this.Boundary.Radial.R(i-1).Outer;
                    R = sparse(1:numel(R), 1:numel(R), exp(1i * R * t));
                    
                    linearMatrix{i,i-1} ...
                        = real(this.Boundary.Radial.S(i).Inner * this.Boundary.Radial.G(i-1).Outer * R * this.Boundary.Radial.D(i-1).Outer);
                end
                
                if i < nAssemblies
                    R = this.Boundary.Radial.R(i).Outer - this.Boundary.Radial.R(i+1).Inner;
                    R = sparse(1:numel(R), 1:numel(R), exp(1i * R * t));
                    
                    linearMatrix{i,i+1} ...
                        = real(this.Boundary.Radial.S(i).Outer * this.Boundary.Radial.G(i+1).Inner * R * this.Boundary.Radial.D(i+1).Inner);
                end
            end
            linearMatrix = cell2mat(linearMatrix);
        end
        
        function conductivityMatrix = C(this, ~, h)
            assembly           = this.Assemblies_;
            nAssemblies        = numel(assembly);
            conductivityMatrix = cell(nAssemblies,nAssemblies);
           
            for i = 1:nAssemblies
                conductivityMatrix{i,i} =   this.Mass.Conductivity(i).El2El / h ...
                                          + this.Mass.Conductivity(i).Re2El / h...
                                          + this.Mass.Conductivity(i).El2Re ...
                                          + this.Mass.Conductivity(i).Re2Re;
            end

            conductivityMatrix = blkdiag(conductivityMatrix{:});
        end
        
    	function exogenousFunction  = f(this, t, h)
            if nargin < 3
                h = 1;
            end
            
            mesh     = this.Mesh_;
            assembly = [mesh.Assembly];
            
            nAssemblies       = numel(assembly);
            exogenousFunction = cell(nAssemblies,1);
            for i = 1:nAssemblies
                materials = [assembly(i).Regions.Material];
                [~,Mx,My] = materials.vectorMr;
                Mx        = Mx(mesh(i).ElementRegions);
                My        = My(mesh(i).ElementRegions);
                curl      = this.Jacobian.MagnetizationCurrent(i);
                
                exogenousFunction{i,1} = (- curl.dIzdMx * Mx - curl.dIzdMy * My);
                
               	sources = assembly(i).Sources;
                if numel(sources) > 0
                    exogenousFunction{i,1} = exogenousFunction{i,1} + (   this.Exogenous.Magnetic(i).Sr2El     ...
                                                                        + this.Exogenous.Magnetic(i).Sr2Re * h ...
                                                                        + this.Exogenous.Magnetic(i).Sr2Cr * h) * sources.f(t);
                end
            end
            exogenousFunction = cell2mat(exogenousFunction);
        end
        
        %% Postprocessing Functions
        function [y, y_t] = doPostProcessing(this, x, x_t)
            %% Recover Full Solution
            nTimes      = numel(x);
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            y           = cell(nAssemblies,nTimes);
            y_t         = cell(nAssemblies,nTimes);
            index       = this.Index;
            for i = 1:nAssemblies
                I = index.Global(i).X;
                for j = 1:nTimes
                    y{i,j}   = ppMatrices(i).Reduced2Full * x{j}(I);
                    y_t{i,j} = ppMatrices(i).Reduced2Full * x_t{j}(I);
                end
            end
        end
    end
end