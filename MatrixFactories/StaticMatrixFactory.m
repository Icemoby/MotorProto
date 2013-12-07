classdef StaticMatrixFactory < MatrixFactory
    methods
        function this = StaticMatrixFactory(varargin)
            warning('MotorProto:Verbose', 'Current implementation assumes anti-periodic symmetry');
            if nargin > 0
                this.Model = varargin{1};
                this       = build(this, 'space');
            end
        end
        
        %% preprocessing functions
    	function this = buildPostProcessingMatrices(this)
            warning('MotorProto:Verbose', 'Current model does not verify that all assembly lengths are specified identically. Generalize length calculation');
            warning('MotorProto:Verbose', 'Current implementation assumes anti-periodic symmetry');
            
            mesh           = this.Mesh_;
            nMesh          = numel(mesh);
            index          = this.Index;
            
            postProcessing = struct.empty(0,nMesh);
            for i = 1:nMesh;
                assembly     = mesh(i).Assembly;
                components   = assembly.Components;
                regions      = assembly.Regions;
                sources      = assembly.Sources;
                
                nUnknowns    = index.Local(i).Unknowns;
                nComponents  = numel(components);
                nSources     = numel(sources);
                
                I            = speye(nUnknowns, nUnknowns);
                
                elements     = mesh(i).Elements;
                elRegion     = mesh(i).ElementRegions;
                elAreas      = mesh(i).ElementAreas;

                nElements = numel(elRegion);
                
                assert(nSources == 0 || nSources == 1, 'MotorProto:MatrixFactory', 'Only one source per assembly is currently supported');
                
                %% Boundary Reduction Matrices
              	postProcessing(i).Full2Reduced = this.applyPeriodicBoundaryConditions(I, i, 'ReduceRow');
                postProcessing(i).Reduced2Full = this.applyPeriodicBoundaryConditions(I, i, 'RestoreRow');
                
                %% Solution to Magnetic Vector Potential Conversion
                postProcessing(i).X2A = I(index.Local(i).A,:);

                %% Source Contribution Matrices
                if nSources > 0
                    nPhases            = sources.Phases.Value;
                    connectionMatrices = sources.ConnectionMatrices;
                    connectionPolarity = sources.ConnectionPolarity;
                    
                    phaseResistance    = zeros(1,nPhases);
                    maxEntries         = 0;
                    
                    %% Calculate Resistances
                    for k = 1:nPhases
                        [nSeries, nParallel]   = size(connectionMatrices{k});
                        maxEntries             = maxEntries + nSeries * nParallel * nPhases;
                        connectionRegions      = components(connectionMatrices{k});
                        connectionGeometry     = [connectionRegions.Geometry];
                        connectionMaterial     = [connectionRegions.Material];
                        connectionConductivity = [connectionMaterial.sigma];
                        connectionArea         = [connectionGeometry.area];
                        connectionResistance   = assembly.Length.Value ./ (connectionConductivity .* connectionArea);
                        connectionResistance   = reshape(connectionResistance, nSeries, nParallel);
                        
                        seriesResistance       = sum(connectionResistance, 1);
                        phaseResistance(k)     = 1 ./ sum(1 ./ seriesResistance, 2);
                    end
                    
                    rows = zeros(1, maxEntries);
                    cols = zeros(1, maxEntries);
                    valJ = zeros(1, maxEntries);
                    valE = zeros(1, maxEntries);
                    F2V  = zeros(nPhases, nPhases);
                    F2I  = zeros(nPhases, nPhases);
                    ind  = 0;
                    
                    %% Calculate Matrix
                    if (sources.ConnectionType == ConnectionTypes.Wye)
                        parallelResistance = 1 ./ sum(1 ./ phaseResistance, 2);
                        for k = 1:nPhases
                            [nSeries, nParallel]   = size(connectionMatrices{k});
                            connectionRegions      = components(connectionMatrices{k});
                            connectionMaterial     = [connectionRegions.Material];
                            connectionConductivity = [connectionMaterial.sigma];
                            
                            if (sources.Type == SourceTypes.VoltageSource) && (assembly.ConnectionType == ConnectionTypes.Wye)
                                %% Wye Connected Voltage Source, Wye Connected Load
                                %   Electric Field, Current Density
                                for l = 1:nPhases
                                    I = (ind + 1):(ind + nSeries * nParallel);
                                    
                                    rows(I) = reshape(connectionMatrices{k}, 1, []);
                                    cols(I) = l;
                                    valE(I) = reshape(connectionPolarity{k}, 1, []);
                                    valE(I) = valE(I) / (assembly.Length.Value * nSeries) * assembly.ModeledFraction;
                                    if k ~= l
                                        valE(I) = valE(I) * (  - parallelResistance / phaseResistance(l));
                                    else
                                        valE(I) = valE(I) * (1 - parallelResistance / phaseResistance(l));
                                    end
                                    
                                    valJ(I) = valE(I) .* connectionConductivity;
                                    
                                    ind = ind + nSeries * nParallel;
                                end

                                %   Voltage
                                F2V(k,k) = 1;
                                F2V(k,:) = F2V(k,:) - parallelResistance ./ phaseResistance;
                                
                                %   Current
                                F2V(k,:) = F2V(k,:) / phaseResistance(k) * assembly.ModeledFraction;
                            elseif (sources.Type == SourceTypes.VoltageSource) && (assembly.ConnectionType == ConnectionTypes.Delta)
                                %% Wye Connected Voltage Source, Delta Connected Load
                                %   Electric Field, Current Density
                                I = (ind + 1):(ind + nSeries * nParallel);
                                
                                rows(I) = reshape(connectionMatrices{k}, 1, []);
                                cols(I) = k;
                                valE(I) = reshape(connectionPolarity{k}, 1, []);
                                valE(I) = valE(I) / (assembly.Length.Value * nSeries) * assembly.ModeledFraction;
                                
                                valJ(I) = valE(I) .* connectionConductivity;
                                
                                ind     = ind + nSeries * nParallel;
                                
                                I       = (ind + 1):(ind + nSeries * nParallel);
                                
                                warning('MotorProto:Verbose', 'This can be generalized to allow all possible m-phase wye-delta connections');
                                l       = mod(k + ceil(nPhases / 2) - 1, nPhases) + 1; %maximum voltage connection
                                
                                rows(I) = reshape(connectionMatrices{k}, 1, []);
                                cols(I) = l;
                                valE(I) = - reshape(connectionPolarity{k}, 1, []);
                                valE(I) = valE(I) / (assembly.Length.Value * nSeries) * assembly.ModeledFraction;
                                
                                valJ(I) = valE(I) .* connectionConductivity;
                                
                                ind     = ind + nSeries * nParallel;
                                
                                %   Voltage
                                F2V(k,k) =  1;
                                F2V(k,l) = -1;
                                
                                %   Current
                                F2I(k,:) = F2V(k,:) / phaseResistance(k) * assembly.ModeledFraction;
                            else
                                error('MotorProto:MatrixFactory', 'No implementation for %s connected %s with a %s connected load',...
                                    char(sources.ConnectionType), char(sources.Type), char(assembly.ConnectionType));
                            end
                        end
                    elseif (sources.ConnectionType == ConnectionTypes.Delta)
                        seriesResistance = sum(phaseResistance, 2);
                        for k = 1:nPhases
                            [nSeries, nParallel]   = size(connectionMatrices{k});
                            connectionRegions      = components(connectionMatrices{k});
                            connectionMaterial     = [connectionRegions.Material];
                            connectionConductivity = [connectionMaterial.sigma];
                            
                            if (sources.Type == SourceTypes.CurrentSource) && (assembly.ConnectionType == ConnectionTypes.Delta)
                                %% Delta Connected Current Source, Delta Connected Load
                                %   Current Density, Electric Field
                                for l = 1:nPhases
                                    I = (ind + 1):(ind + nSeries * nParallel);
                                    
                                    rows(I) = reshape(connectionMatrices{k}, 1, []);
                                    cols(I) = l;
                                    valJ(I) = reshape(connectionPolarity{k}, 1, []);
                                    valJ(I) = valJ(I) .* (connectionConductivity / assembly.Length.Value / nSeries);
                                    if k ~= l
                                        valJ(I) = valJ(I) * (  - phaseResistance(l) / seriesResistance) * phaseResistance(k);
                                    else
                                        valJ(I) = valJ(I) * (1 - phaseResistance(l) / seriesResistance) * phaseResistance(k);
                                    end
                                    
                                    valE(I) = valJ(I) ./ connectionConductivity;
                                    ind = ind + nSeries * nParallel;
                                end
                                
                                %	Current
                                F2I(k,k) = 1;
                                F2I(k,:) = F2I(k,:) - phaseResistance / seriesResistance;
                                
                                %   Voltage
                                F2V(k,:) = F2I(k,:) * phaseResistance(k) / assembly.ModeledFraction;
                            elseif (sources.Type == SourceTypes.CurrentSource) && (assembly.ConnectionType == ConnectionTypes.Wye)
                                %% Delta Connected Current Source, Wye Connected Load
                                %   Current Density, Electric Field
                                I = (ind + 1):(ind + nSeries * nParallel);
                                
                                rows(I) = reshape(connectionMatrices{k}, 1, []);
                                cols(I) = k;
                                valJ(I) = reshape(connectionPolarity{k}, 1, []);
                                valJ(I) = valJ(I) .* (connectionConductivity / assembly.Length.Value / nSeries) *  phaseResistance(k);
                              	valE(I) = valJ(I) ./ connectionConductivity;
                                
                                ind     = ind + nSeries * nParallel;
                                
                                I       = (ind + 1):(ind + nSeries * nParallel);
                                
                                warning('MotorProto:Verbose', 'This can be generalized to allow all possible m-phase delta-wye connections');
                                l       = mod(k + ceil(nPhases / 2) - 1, nPhases) + 1; %maximum current connection
                                
                                rows(I) = reshape(connectionMatrices{k}, 1, []);
                                cols(I) = l;
                                valJ(I) = - reshape(connectionPolarity{k}, 1, []);
                                valJ(I) = valJ(I) .* (connectionConductivity / assembly.Length.Value / nSeries) * phaseResistance(l);
                                valE(I) = valJ(I) ./ connectionConductivity;
                                
                                ind     = ind + nSeries * nParallel;
                                
                                %   Current
                                F2I(k,k) = 1;
                                F2I(k,l) = -1;
                                
                                %   Voltage
                                F2V(k,:) = phaseResistance(k) * F2I(k,:) / assembly.ModeledFraction;
                            else
                                error('MotorProto:MatrixFactory', 'No implementation for %s connected %s with a %s connected load',...
                                    char(sources.ConnectionType), char(sources.Type), char(assembly.ConnectionType));
                            end
                        end
                    else
                        error('MotorProto:MatrixFactory','No implementation for %s connected sources',char(sources.ConnectionType));
                    end
                    
                    rows(ind+1:end) = [];
                    cols(ind+1:end) = [];
                    valJ(ind+1:end) = [];
                    valE(ind+1:end) = [];
                    
                    %% Source to Component Current Density
                    Sr2J = sparse(rows, cols, valJ, nComponents, nPhases);
                    Sr2E = sparse(rows, cols, valE, nComponents, nPhases);
                    
                    %% Component to Element Conversion
                    ind = 0;
                    nEntries = length(mesh(i).Elements);
                    rows     = zeros(1, nEntries);
                    cols     = zeros(1, nEntries);
                    valJ     = zeros(1, nEntries);
                    for k = 1:nComponents
                        kRegion = find(components(k) == regions);
                        if ~isempty(kRegion)
                            iElement = find(mesh(i).ElementRegions == kRegion);
                            
                            I        = (ind+1):(ind+numel(iElement));
                            
                            rows(I)  = iElement;
                            cols(I)  = k;
                            valJ(I)  = 1;
                            
                            ind      = ind + numel(iElement);
                        end
                    end
                    
                    C2El = sparse(rows, cols, valJ, length(mesh(i).Elements), nComponents);
                    
                    %% Source to Element
                    F2J = C2El * Sr2J;
                    F2E = C2El * Sr2E;

                    postProcessing(i).F2J = F2J;
                    postProcessing(i).F2E = F2E;
                    postProcessing(i).F2V = F2V;
                    postProcessing(i).F2I = F2I;
                    
                else
                    postProcessing(i).F2J = sparse(length(mesh(i).Elements), 0);
                    postProcessing(i).F2E = sparse(length(mesh(i).Elements), 0);
                    postProcessing(i).F2V = sparse(length(mesh(i).Elements), 0);
                    postProcessing(i).F2I = sparse(length(mesh(i).Elements), 0);
                end
                
                %% Solution Contribution Matrices
                %	Electric Field
                postProcessing(i).X2E          = sparse(nElements,nUnknowns);
                [postProcessing(i).X_t2E{1:3}] = deal(sparse(nElements,nUnknowns));  
                
                %   Current Density
                postProcessing(i).X2J          = sparse(nElements,nUnknowns);
                [postProcessing(i).X_t2J{1:3}] = deal(sparse(nElements,nUnknowns));
                
                if nSources > 0
                    nPhases = sources.Phases.Value;
                    %   Voltage
                    postProcessing(i).X2V   = sparse(nPhases, nUnknowns);
                    
                    %   Current
                    postProcessing(i).X2I   = sparse(nPhases, nUnknowns);
                    postProcessing(i).X_t2I = sparse(nPhases, nUnknowns);
                    
                    %   Flux Linkage
                    rows = zeros(1,0);
                    cols = zeros(1,0);
                    vals = zeros(1,0);
                    
                    connectionMatrices = sources.ConnectionMatrices;
                    connectionPolarity = sources.ConnectionPolarity;
                    ind2 = 0;
                    for j = 1:nPhases
                        connectionComponents = components(connectionMatrices{j});
                        [nSeries, nParallel] = size(connectionMatrices{j});
                        
                        if ~all(size(connectionComponents) == size(connectionMatrices{j}))
                            connectionComponents = connectionComponents.';
                        end
                        
                        for k = 1:nSeries
                            ind1      = ind2;
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
                                
                                rows(ind2+1:ind2+3*N) = j;
                                cols(ind2+1:ind2+3*N) = nodes;
                                vals(ind2+1:ind2+3*N) = connectionPolarity{j}(k,l) * areas * assembly.Length.Value / (3 * assembly.ModeledFraction);
                                
                                ind2 = ind2 + 3*N;
                            end
                            vals(ind1+1:ind2) = vals(ind1+1:ind2) / totalArea;
                        end
                    end
                    
                    postProcessing(i).X2Lambda = sparse(rows, cols, vals, nPhases, nUnknowns);
                    postProcessing(i).X_t2V    = sparse(rows, cols, vals, nPhases, nUnknowns);
                else
                    %   Voltage
                    postProcessing(i).X2V   = sparse(0, nUnknowns);
                    
                    %   Current
                    postProcessing(i).X2I   = sparse(0, nUnknowns);
                    postProcessing(i).X_t2I = sparse(0, nUnknowns);
                    
                    %   Flux Linkage
                    postProcessing(i).X2Lambda = sparse(0, nUnknowns);
                    postProcessing(i).X_t2V    = sparse(0, nUnknowns);
                    
                end
            end
            this.PostProcessing = postProcessing;
        end
        
        %% Solver functions
        function                           linearMatrix = K(this, t, h)
            warning('MotorProto:Verbose','Change the form of this function to be more like the individual component creation methods');
            %% create stiffness matrix
            %  periodic boundary conditions are included in the curl matrices
            assembly     = this.Assemblies_;
            nAssemblies  = numel(assembly);
            linearMatrix = cell(nAssemblies, nAssemblies);
            
            for i = 1:nAssemblies
                linearMatrix{i,i} = this.Stiffness.Reluctivity(i).El2El;
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
        
        function [nonlinearJacobian, nonlinearFunction] = G(this, t, x)
            mesh     = this.Mesh_;
            assembly = [mesh.Assembly];
            index    = this.Index;
            curlE2N  = this.Jacobian.MagnetizationCurrent;
            curlN2E  = this.Jacobian.FluxDensity;
            
            nMesh     = numel(mesh);
            nUnknowns = index.Global(end).Unknowns;
            nRows     = zeros(nMesh,1);
            for i = 1:nMesh
                nRows(i) = index.Local(i).Unknowns;
            end
            
            nonlinearFunction = mat2cell(sparse(nUnknowns,1), nRows, 1);
            nonlinearJacobian = mat2cell(sparse(nUnknowns,nUnknowns), nRows, nRows);

            for i = 1:nMesh
                materials = [assembly(i).Regions.Material];
                I         = index.Global(i).X;
                
                Bx     = curlN2E(i).dBxdXz * x(I);
                By     = curlN2E(i).dBydXz * x(I);
                sizeB  = size(Bx);
                Mx     = zeros(sizeB);
                My     = zeros(sizeB);
                dMxdBx = zeros(sizeB);
                dMydBy = zeros(sizeB);
                dMydBx = zeros(sizeB);
                dMxdBy = zeros(sizeB);

                for j = 1:numel(materials)
                    if ~materials(j).Linear
                        J = (mesh(i).ElementRegions == j);
                        [Mx(J) ,My(J), dMxdBx(J), dMydBy(J), dMydBx(J), dMxdBy(J)] = materials(j).vectorM(Bx(J), By(J));
%                         %% Linear Model Debugging Code
%                         Mx(J,:) = Bx(J,:)*(1-1/1000)/mu_o;
%                         My(J,:) = By(J,:)*(1-1/1000)/mu_o;
%                         dMxdBx(J,:) = (1-1/1000)/mu_o;
%                         dMydBy(J,:) = (1-1/1000)/mu_o;
%                         dMydBx(J,:) = 0;
%                         dMxdBy(J,:) = 0;
                    end
                end
                
                I = 1:sizeB(1);
                nonlinearJacobian{i,i} ...
                    = + curlE2N(i).dIzdMx*( sparse(I,I,dMxdBx) * curlN2E(i).dBxdXz + sparse(I,I,dMxdBy) * curlN2E(i).dBydXz)...
                      + curlE2N(i).dIzdMy*( sparse(I,I,dMydBx) * curlN2E(i).dBxdXz + sparse(I,I,dMydBy) * curlN2E(i).dBydXz);

                nonlinearFunction{i} = + curlE2N(i).dIzdMx * Mx  + curlE2N(i).dIzdMy * My;
            end
            
            nonlinearJacobian = cell2mat(nonlinearJacobian);
            nonlinearFunction = cell2mat(nonlinearFunction);
        end
        
        function                     conductivityMatrix = C(this, t, h)
        	conductivityMatrix = this.Conductivity.El2El;
        end
        
        function                      exogenousFunction = f(this, t, h)
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
                exo       = this.Exogenous.Magnetic(i);
                
                exogenousFunction{i,1} = - curl.dIzdMx * Mx  - curl.dIzdMy * My;
                
               	sources = assembly(i).Sources;
                if numel(sources) > 0
                    exogenousFunction{i, 1} = exogenousFunction{i, 1} + exo.Sr2El * sources.f(t);
                end
            end
            exogenousFunction = cell2mat(exogenousFunction);
        end
        
        %% Postprocessing functions
        function [y, y_t] = doPostProcessing(this, x, ~)
            if ~iscell(x)
                [n,m] = size(x);
                x     = mat2cell(x, n, ones(1,m));
            end
            
            nTimes = numel(x);
            
            %% Estimate Derivatives
            if nTimes > 1
                x_t    = x(:, 1:(nTimes - 1));
                x_t    = cell2mat(x_t);
                fe     = [this.Assemblies_.ElectricalFrequency];
                fe     = [fe.Value];

                sameFrequency = (abs(fe - mean(fe)) < sqrt(eps));
                assert(all(sameFrequency), 'MotorProto:StaticMatrixFactory', 'All electrical frequencies must be the same');
                fe           = mean(fe);
                hMax         = (nTimes - 1) / 2 - 1;
                h            = [1:(hMax+1), (hMax+3):(2*hMax + 2)];
                nyqEl        = hMax + 2;
                I            = 1:(2*hMax + 1);
                omega        = 1i * 2 * pi * fe * sparse(I, I, [0:1:hMax -hMax:1:-1]);
                x_t          = fft(x_t, [], 2);
                x_t(:,h)     = x_t(:,h) * omega;
                x_t(:,nyqEl) = 0;
                x_t          = ifft(x_t,[],2);
                x_t          = mat2cell(x_t, length(x_t(:,1)), ones(1,nTimes - 1));
                x_t(:,end+1) = x_t(:,1);
            else
                x_t    = x;
                x_t{1} = x_t{1} * 0;
            end
            
            %% Recover Full Solution
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            y           = cell(nAssemblies,nTimes);
            y_t         = cell(nAssemblies,nTimes);
            index       = this.Index;
            for i = 1:nAssemblies
                I = index.Global(i).X;
                for j = 1:nTimes
                    y  {i, j} = ppMatrices(i).Reduced2Full * x  {j}(I);
                    y_t{i, j} = ppMatrices(i).Reduced2Full * x_t{j}(I);
                end
            end
        end
    end
end