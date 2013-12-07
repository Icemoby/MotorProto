classdef HarmonicMatrixFactory < DynamicMatrixFactory
	methods
        function this = HarmonicMatrixFactory(varargin)
            warning('MotorProto:Verbose',...
                    'Present implementation assumes anti-periodic symmetry');
            if nargin > 0
                this.Model = varargin{1};
                this = build(this,'spacetime');
%                 this = build(this,'space');
            end
        end

        %% Preprocessing Functions
        function [structure, newFields, pbcType] = buildLocalIndexVectors(this,structure, iMesh)
            [structure, newFields, pbcType] = buildLocalIndexVectors@DynamicMatrixFactory(this, structure, iMesh);
        end
        
        %% Boundary Matrices
        function [structure, newFields, pbcType] = buildRadialBoundaryMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildRadialBoundaryMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            for i = 1:numel(newFields)
                subFields = fields(structure.(newFields{i}));
                for j = 1:numel(subFields)
                    if ~strcmpi(pbcType{i}, 'none')
                        structure.(newFields{i})(iMesh).(subFields{j}) ...
                            = {     structure.(newFields{i})(iMesh).(subFields{j}), 0 * structure.(newFields{i})(iMesh).(subFields{j});
                                0 * structure.(newFields{i})(iMesh).(subFields{j}),     structure.(newFields{i})(iMesh).(subFields{j})};
                    end
                end
            end
        end
        
        %% Stiffness Matrices
        function [structure, newFields, pbcType] = buildReluctivityMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildReluctivityMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            reluctivity = structure.Reluctivity;
            
            reluctivity(iMesh).Full =   reluctivity(iMesh).El2El + reluctivity(iMesh).Re2El ...
                                      + reluctivity(iMesh).El2Re + reluctivity(iMesh).Re2Re ...
                                      + reluctivity(iMesh).Re2Cr + reluctivity(iMesh).Cr2Re ...
                                      + reluctivity(iMesh).Cr2Cr;
            
            flds = fields(reluctivity);
            for i = 1:numel(flds)
                reluctivity(iMesh).(flds{i}) = {    reluctivity(iMesh).(flds{i}), 0 * reluctivity(iMesh).(flds{i});
                                                0 * reluctivity(iMesh).(flds{i}),     reluctivity(iMesh).(flds{i})};
            end
                                  
            structure.Reluctivity = reluctivity;
        end
        
        %% Mass Matrices
        function [structure, newFields, pbcType] = buildConductivityMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildConductivityMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            conductivity = structure.Conductivity;
            
            conductivity(iMesh).Full =   conductivity(iMesh).El2El + conductivity(iMesh).Re2El ...
                                       + conductivity(iMesh).El2Re + conductivity(iMesh).Re2Re;
                            
            flds = fields(conductivity);
            for i = 1:numel(flds)
                conductivity(iMesh).(flds{i}) = {0 * conductivity(iMesh).(flds{i}), -1 * conductivity(iMesh).(flds{i});
                                                     conductivity(iMesh).(flds{i}),  0 * conductivity(iMesh).(flds{i})};
            end
                
                                   
            structure.Conductivity = conductivity;
        end
        
        %% Jacobian Matrices
      	function [structure, newFields, pbcType] = buildMagnetizationJacobianMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildMagnetizationJacobianMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            magnetizationCurrent = structure.MagnetizationCurrent;
            
            flds = fields(magnetizationCurrent);
            for i = 1:numel(flds)
                magnetizationCurrent(iMesh).(flds{i}) = {    magnetizationCurrent(iMesh).(flds{i}), 0 * magnetizationCurrent(iMesh).(flds{i});
                                                         0 * magnetizationCurrent(iMesh).(flds{i}),     magnetizationCurrent(iMesh).(flds{i})};
            end
                                  
            structure.MagnetizationCurrent = magnetizationCurrent;
        end
        
        function [structure, newFields, pbcType] = buildFluxDensityJacobianMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildFluxDensityJacobianMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            fluxDensity = structure.FluxDensity;
            
            flds = fields(fluxDensity);
            for i = 1:numel(flds)
                fluxDensity(iMesh).(flds{i}) = {    fluxDensity(iMesh).(flds{i}), 0 * fluxDensity(iMesh).(flds{i});
                                                0 * fluxDensity(iMesh).(flds{i}),     fluxDensity(iMesh).(flds{i})};
            end
                                  
            structure.FluxDensity = fluxDensity;
        end
        
        %% Exogenous Matrices
        function [structure, newFields, pbcType] = buildMagneticInputMatrices(this, structure, iMesh)
            [structure, newFields, pbcType] = buildMagneticInputMatrices@DynamicMatrixFactory(this, structure, iMesh);
            
            magnetic = structure.Magnetic;
            
            magnetic(iMesh).Full = magnetic(iMesh).Sr2El + magnetic(iMesh).Sr2Re + magnetic(iMesh).Sr2Cr;
            
            flds = fields(magnetic);
            for i = 1:numel(flds)
                magnetic(iMesh).(flds{i}) = {    magnetic(iMesh).(flds{i}), 0 * magnetic(iMesh).(flds{i});
                                             0 * magnetic(iMesh).(flds{i}),     magnetic(iMesh).(flds{i})};
            end
            
            structure.Magnetic = magnetic;
        end
        
        %% Boundary Conditions
        function mats = applyPeriodicBoundaryConditions(this, mats, iMesh, action)
            if ~strcmpi(action,'none')
                mesh  = this.Mesh_(iMesh);
                I     = mesh.PeriodicBoundaryNodes(1, :);
                J     = mesh.PeriodicBoundaryNodes(2, :);
                N     = numel(mesh.X);
                K     = setdiff(1:N, [I, J]);
%                 C     = -eye(2);
                coef  = mesh.Assembly.SolutionSpaceTimeCoefficients;
                c     = [cos(2 * pi / coef(3)) -sin(2 * pi / coef(3));
                         sin(2 * pi / coef(3))  cos(2 * pi / coef(3))].';
                     
                if isstruct(mats)
                    field = fields(mats);
                    for i = 1:numel(field)                   
                        %% Duplicate for the number of additional space/time symmetries
                        mats(iMesh).(field{i}) = repmat({mats(iMesh).(field{i})}, 1, coef(2));
                        
                        switch lower(action)
                            case {'row', 'rows'}
                                for j = 1:coef(2)
                                	C = coef(1) * c^(j-1);
                                    M = mats(iMesh).(field{i}){j};
                                    N = M;
                                    
                                    %% Real Part (cos)
                                  	M{1,1}(I,:) = M{1,1}(I,:) + C(1,1) * N{1,1}(J,:);
                                    M{2,1}(I,:) = M{2,1}(I,:) + C(1,1) * N{2,1}(J,:);
                                    M{1,2}(I,:) = M{1,2}(I,:) + C(2,1) * N{1,1}(J,:);
                                    M{2,2}(I,:) = M{2,2}(I,:) + C(2,1) * N{2,1}(J,:);
                                    
                                    %% Imag Part (sin)
                                    M{2,2}(I,:) = M{2,2}(I,:) + C(2,2) * N{2,2}(J,:);
                                    M{1,2}(I,:) = M{1,2}(I,:) + C(2,2) * N{1,2}(J,:);
                                    M{2,1}(I,:) = M{2,1}(I,:) + C(1,2) * N{2,2}(J,:);
                                    M{1,1}(I,:) = M{1,1}(I,:) + C(1,2) * N{1,2}(J,:);
                                    
                                    %% Delete Rows
                                    M{1,1}(J,:) = [];
                                    M{1,2}(J,:) = [];
                                    M{2,1}(J,:) = [];
                                    M{2,2}(J,:) = [];
                                    
                                    %% Update Structure
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            case {'column', 'columns'}
                                for j = 1:coef(2)
                                	C = coef(1) * c^(j-1);
                                    M = mats(iMesh).(field{i}){j};
                                    N = M;
                                    
                                    %% Real Part (cos)
                                  	M{1,1}(:,I) = M{1,1}(:,I) + C(1,1) * N{1,1}(:,J);
                                    M{2,1}(:,I) = M{2,1}(:,I) + C(1,1) * N{2,1}(:,J);
                                    M{1,2}(:,I) = M{1,2}(:,I) + C(1,2) * N{1,1}(:,J);
                                    M{2,2}(:,I) = M{2,2}(:,I) + C(1,2) * N{2,1}(:,J);
                                    
                                    %% Imag Part (sin)
                                    M{2,2}(:,I) = M{2,2}(:,I) + C(2,2) * N{2,2}(:,J);
                                    M{1,2}(:,I) = M{1,2}(:,I) + C(2,2) * N{1,2}(:,J);
                                    M{2,1}(:,I) = M{2,1}(:,I) + C(2,1) * N{2,2}(:,J);
                                    M{1,1}(:,I) = M{1,1}(:,I) + C(2,1) * N{1,2}(:,J);
                                    
                                    %% Delete Columns
                                    M{1,1}(:,J) = [];
                                    M{1,2}(:,J) = [];
                                    M{2,1}(:,J) = [];
                                    M{2,2}(:,J) = [];
                                    
                                    %% Update Structure
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            case 'both'
                                for j = 1:coef(2)
                                    M = mats(iMesh).(field{i}){j};
                                    N = M;
                                    
                                    %% Adjust Diagonal Components
                                    M{1,1}(I,I) = M{1,1}(I,I) + N{1,1}(J,J);
                                    M{1,2}(I,I) = M{1,2}(I,I) + N{1,2}(J,J);
                                    M{2,1}(I,I) = M{2,1}(I,I) + N{2,1}(J,J);
                                    M{2,2}(I,I) = M{2,2}(I,I) + N{2,2}(J,J);
                                    
                                    %% Apply Boundary Conditions
                                    %% Rows
                                	C = coef(1) * c^(j-1);
                                    %% Real Part (cos)
                                  	M{1,1}(I,K) = M{1,1}(I,K) + C(1,1) * N{1,1}(J,K);
                                    M{2,1}(I,K) = M{2,1}(I,K) + C(1,1) * N{2,1}(J,K);
                                    M{1,2}(I,K) = M{1,2}(I,K) + C(2,1) * N{1,1}(J,K);
                                    M{2,2}(I,K) = M{2,2}(I,K) + C(2,1) * N{2,1}(J,K);
                                    
                                    %% Imag Part (sin)
                                    M{2,2}(I,K) = M{2,2}(I,K) + C(2,2) * N{2,2}(J,K);
                                    M{1,2}(I,K) = M{1,2}(I,K) + C(2,2) * N{1,2}(J,K);
                                    M{2,1}(I,K) = M{2,1}(I,K) + C(1,2) * N{2,2}(J,K);
                                    M{1,1}(I,K) = M{1,1}(I,K) + C(1,2) * N{1,2}(J,K);
                                    
                                    %% Columns
                                	C = coef(1) * c^(j-1);
                                    %% Real Part (cos)
                                  	M{1,1}(K,I) = M{1,1}(K,I) + C(1,1) * N{1,1}(K,J);
                                    M{2,1}(K,I) = M{2,1}(K,I) + C(1,1) * N{2,1}(K,J);
                                    M{1,2}(K,I) = M{1,2}(K,I) + C(1,2) * N{1,1}(K,J);
                                    M{2,2}(K,I) = M{2,2}(K,I) + C(1,2) * N{2,1}(K,J);
                                    
                                    %% Imag Part (sin)
                                    M{2,2}(K,I) = M{2,2}(K,I) + C(2,2) * N{2,2}(K,J);
                                    M{1,2}(K,I) = M{1,2}(K,I) + C(2,2) * N{1,2}(K,J);
                                    M{2,1}(K,I) = M{2,1}(K,I) + C(2,1) * N{2,2}(K,J);
                                    M{1,1}(K,I) = M{1,1}(K,I) + C(2,1) * N{1,2}(K,J);
                                    
                                    %% Delete Rows
                                    M{1,1}(J,:) = [];
                                    M{1,2}(J,:) = [];
                                    M{2,1}(J,:) = [];
                                    M{2,2}(J,:) = [];
                                    
                                    %% Delete Columns
                                    M{1,1}(:,J) = [];
                                    M{1,2}(:,J) = [];
                                    M{2,1}(:,J) = [];
                                    M{2,2}(:,J) = [];
                                    
                                    %% Update Structure
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            case 'reducerow'
                                for j = 1:coef(2)
                                	M = mats(iMesh).(field{i}){j};
                                    
                                    M{1,1}(J,:) = [];
                                    M{1,2}(J,:) = [];
                                    M{2,1}(J,:) = [];
                                    M{2,2}(J,:) = [];
                                    
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            case 'restorerow'
                                for j = 1:coef(2)
                                	C = coef(1) * c^(j-1);
                                    M = mats(iMesh).(field{i}){j};
                                    N = M;
                                    
                                    M{1,1}(J,:) = C(1,1) * N{1,1}(I,:);
                                    M{2,1}(J,:) = C(2,1) * N{1,1}(I,:);
                                    M{1,2}(J,:) = C(1,2) * N{2,2}(I,:);
                                    M{2,2}(J,:) = C(2,2) * N{2,2}(I,:);
                                    
                                    M{1,1}(:,J) = [];
                                    M{1,2}(:,J) = [];
                                    M{2,1}(:,J) = [];
                                    M{2,2}(:,J) = [];
                                    
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            case 'boundaryoperation'
                                for j = 1:coef(2)
                                    nSections = 1 / mesh.Assembly.ModeledFraction;
                                    harmonics = this.getRadialBoundaryHarmonics(iMesh);
                                    switch lower(field{i})
                                        case 'inner'
                                            harmonics = harmonics{1};
                                        case 'outer'
                                            harmonics = harmonics{2};
                                    end
                                    nHarmonics = numel(harmonics);
                                    
                                    %% Account for each unmodeled region using space/time symmetry
                                    M = mats(iMesh).(field{i}){j};
                                   	N = M;
                                    C = coef(1) * c^(j-1);
                                    
                                    %   Real Part (cos)
                                  	M{1,1}(:,I) = M{1,1}(:,I) + C(1,1) * N{1,1}(:,J);
                                    M{2,1}(:,I) = M{2,1}(:,I) + C(1,1) * N{2,1}(:,J);
                                    M{1,2}(:,I) = M{1,2}(:,I) + C(1,2) * N{1,1}(:,J);
                                    M{2,2}(:,I) = M{2,2}(:,I) + C(1,2) * N{2,1}(:,J);
                                    
                                    %   Imag Part (sin)
                                    M{2,2}(:,I) = M{2,2}(:,I) + C(2,2) * N{2,2}(:,J);
                                    M{1,2}(:,I) = M{1,2}(:,I) + C(2,2) * N{1,2}(:,J);
                                    M{2,1}(:,I) = M{2,1}(:,I) + C(2,1) * N{2,2}(:,J);
                                    M{1,1}(:,I) = M{1,1}(:,I) + C(2,1) * N{1,2}(:,J);
                                    
                                    %   Remove Nodes
                                    M{1,1}(:,J) = [];
                                    M{1,2}(:,J) = [];
                                    M{2,1}(:,J) = [];
                                    M{2,2}(:,J) = [];
                                    
                                    R = sparse(1:nHarmonics, 1:nHarmonics, exp(-1i * 2 * pi * harmonics / nSections));
                                    N = M;
                                    for k = 1:(nSections-1)
                                        C = (coef(1) * c^(j-1))^k;
                                        
                                        M{1,1} = M{1,1} + R^k * C(1,1) * N{1,1};
                                        M{1,1} = M{1,1} + R^k * C(1,2) * N{2,1};
                                        
                                        M{1,2} = M{1,2} + R^k * C(1,1) * N{1,2};
                                        M{1,2} = M{1,2} + R^k * C(1,2) * N{2,2};
                                        
                                        M{2,1} = M{2,1} + R^k * C(2,1) * N{1,1};
                                        M{2,1} = M{2,1} + R^k * C(2,2) * N{2,1};
                                        
                                        M{2,2} = M{2,2} + R^k * C(2,1) * N{1,2};
                                        M{2,2} = M{2,2} + R^k * C(2,2) * N{2,2};
                                    end
                                    
                                    %% Update Structure
                                    mats(iMesh).(field{i}){j} = M;
                                end
                            otherwise
                                error('MotorProto:MatrixFactory', 'Unknown argument "%s" for method applyPeriodicBoundaryConditions', action);
                        end
                    end
                else
                    switch lower(action)
                        case {'row', 'rows'}
                            error('No Implemenetaiton');
                        case {'column', 'columns'}
                            error('No Implemenetaiton');
                        case 'both'
                            error('No Implemenetaiton');
                        case 'reducerow'
                            mats(J,:) = [];
                        case 'restorerow'
                            mats = {    mats, 0 * mats;
                                    0 * mats,     mats};
                            mats = repmat({mats}, 1, coef(2));
                            for j = 1:coef(2)
                                C = coef(1) * c^(j-1);
                                M = mats{j};
                                N = M;
                                
                                M{1,1}(J,:) = C(1,1) * N{1,1}(I,:);
                                M{2,1}(J,:) = C(2,1) * N{1,1}(I,:);
                                M{1,2}(J,:) = C(1,2) * N{2,2}(I,:);
                                M{2,2}(J,:) = C(2,2) * N{2,2}(I,:);
                                
                                M{1,1}(:,J) = [];
                                M{1,2}(:,J) = [];
                                M{2,1}(:,J) = [];
                                M{2,2}(:,J) = [];
                                
                                mats{j} = M;
                            end
                        otherwise
                            error('MotorProto:MatrixFactory', 'Unknown argument "%s" for method applyPeriodicBoundaryConditions', action);
                    end
                end
            end
        end
        
        function this = applyTangentialBoundaryConditions(this)
           	mesh  = this.Mesh_;
            index = this.Index;
            nMesh = numel(mesh);
            for i = 1:nMesh
%                 coef = mesh(i).Assembly.SolutionSpaceTimeCoefficients;
                
%                 n = numel(mesh(i).X);
%                 I = mesh(i).PeriodicBoundaryNodes(1,:);
                J = mesh(i).PeriodicBoundaryNodes(2,:);
%                 K = setdiff(1:n,[I,J]);
                
                
                %% Stiffness
%                 this.Stiffness.Reluctivity(i).Full ...
%                     = repmat({this.Stiffness.Reluctivity(i).Full},2,2);
                    
%                	this.Stiffness.Reluctivity(i).Full ...
%                 	= repmat({this.Stiffness.Reluctivity(i).Full},1,coef(2));
                
                %% Mass
%                 this.Mass.Conductivity(i).Full...
%                     = repmat({this.Mass.Conductivity(i).Full},2,2);
                
%                 this.Mass.Conductivity(i).Full ...
%                  	= repmat({this.Mass.Conductivity(i).Full},1,coef(2));  
                
                %% Exogenous
%                 this.Exogenous.Magnetic(i).Full ...
%                 	= repmat({this.Exogenous.Magnetic(i).Full},2,2);  
                
%                 this.Exogenous.Magnetic(i).Full ...
%                 	= repmat({this.Exogenous.Magnetic(i).Full},1,coef(2));  
                         
                %% Jacobian
%                 this.Jacobian.Magnetization(i).dIzdMx ...
%                     = repmat({this.Jacobian.Magnetization(i).dIzdMx},2,2);
%                 
%                 this.Jacobian.Magnetization(i).dIzdMx ...
%                     = repmat({this.Jacobian.Magnetization(i).dIzdMx},1,coef(2));
%                 
%                 this.Jacobian.Magnetization(i).dIzdMy ...
%                     = repmat({this.Jacobian.Magnetization(i).dIzdMy},2,2); 
%                 
%                 this.Jacobian.Magnetization(i).dIzdMy ...
%                     = repmat({this.Jacobian.Magnetization(i).dIzdMy},1,coef(2));
%                 
%                 this.Jacobian.Magnetization(i).dBxdXz ...
%                     = repmat({this.Jacobian.Magnetization(i).dBxdXz},2,2);
%                 
%                 this.Jacobian.Magnetization(i).dBxdXz ...
%                     = repmat({this.Jacobian.Magnetization(i).dBxdXz},1,coef(2));
%                 
%               	this.Jacobian.Magnetization(i).dBydXz ...
%                     = repmat({this.Jacobian.Magnetization(i).dBydXz},2,2);
%                 
%                 this.Jacobian.Magnetization(i).dBydXz ...
%                     = repmat({this.Jacobian.Magnetization(i).dBydXz},1,coef(2));
                
                %% PostProcessing
%                 this.PostProcessing(i).Reduced2Full ...
%                     = repmat({this.PostProcessing(i).Reduced2Full},2,2);
%                 
%                 this.PostProcessing(i).Reduced2Full ...
%                     = repmat({this.PostProcessing(i).Reduced2Full},1,coef(2));
%                 
%                 for j = 1:2
%                  	this.Boundary.Radial(i).Radius(j).D = repmat({this.Boundary.Radial(i).Radius(j).D},2,2);
%                     
%                     this.Boundary.Radial(i).Radius(j).S = repmat({this.Boundary.Radial(i).Radius(j).S},2,2);
%                     
%                     this.Boundary.Radial(i).Radius(j).D = repmat({this.Boundary.Radial(i).Radius(j).D},1,coef(2));
%                     
%                     this.Boundary.Radial(i).Radius(j).S = repmat({this.Boundary.Radial(i).Radius(j).S},1,coef(2));
%                     
%                     this.Boundary.Radial(i).Radius(j).P = repmat({this.Boundary.Radial(i).Radius(j).P},1,coef(2));
%                 end
                
%                 for j = 1:coef(2)
                    %% x(J) = C' * x(I), x(I) = C * x(J)
                    %C = coef(1) * exp(1i*2*pi*(j-1)/coef(2));
%                     C = coef(1) * [  cos(2*pi*(j-1)/coef(2)),...
%                                     -sin(2*pi*(j-1)/coef(2));...
%                                      sin(2*pi*(j-1)/coef(2)),...
%                                      cos(2*pi*(j-1)/coef(2))];
%                     C = -eye(2);
                    %% Stiffness
%                     this.Stiffness.Reluctivity(i).Full{j}{1,1}(I,I) ...
%                             =   this.Stiffness.Reluctivity(i).Full{j}{1,1}(I,I) ...
%                               + this.Stiffness.Reluctivity(i).Full{j}{1,1}(J,J);
%                           
%                           
%                     this.Stiffness.Reluctivity(i).Full{j}{2,2}(I,I) ...
%                             =   this.Stiffness.Reluctivity(i).Full{j}{2,2}(I,I) ...
%                               + this.Stiffness.Reluctivity(i).Full{j}{2,2}(J,J);
% 
%                     this.Stiffness.Reluctivity(i).Full{j}{1,1}(I,K) ...
%                     	=            this.Stiffness.Reluctivity(i).Full{j}{1,1}(I,K) ...
%                           + C(1,1) * this.Stiffness.Reluctivity(i).Full{j}{1,1}(J,K);
%                     
%                     this.Stiffness.Reluctivity(i).Full{j}{1,2}(:) = 0;
%                     this.Stiffness.Reluctivity(i).Full{j}{1,2}(I,K) ...
%                     	= C(1,2) * this.Stiffness.Reluctivity(i).Full{j}{1,1}(J,K);
%                     
%                     this.Stiffness.Reluctivity(i).Full{j}{2,1}(:) = 0;
%                   	this.Stiffness.Reluctivity(i).Full{j}{2,1}(I,K) ...
%                     	= C(2,1) * this.Stiffness.Reluctivity(i).Full{j}{2,2}(J,K);
%                     
%                    	this.Stiffness.Reluctivity(i).Full{j}{2,2}(I,K) ...
%                     	=            this.Stiffness.Reluctivity(i).Full{j}{2,2}(I,K) ...
%                           + C(2,2) * this.Stiffness.Reluctivity(i).Full{j}{2,2}(J,K);
%                     
%                     
%                  	this.Stiffness.Reluctivity(i).Full{j}{1,1}(K,I) ...
%                     	=            this.Stiffness.Reluctivity(i).Full{j}{1,1}(K,I) ...
%                           + C(1,1) * this.Stiffness.Reluctivity(i).Full{j}{1,1}(K,J);
%                       
%                     this.Stiffness.Reluctivity(i).Full{j}{1,2}(K,I) ...
%                     	= C(2,1) * this.Stiffness.Reluctivity(i).Full{j}{1,1}(K,J);
%                     
%                   	this.Stiffness.Reluctivity(i).Full{j}{2,1}(K,I) ...
%                     	= C(1,2) * this.Stiffness.Reluctivity(i).Full{j}{2,2}(K,J);
%                     
%                    	this.Stiffness.Reluctivity(i).Full{j}{2,2}(K,I) ...
%                     	=            this.Stiffness.Reluctivity(i).Full{j}{2,2}(K,I) ...
%                           + C(2,2) * this.Stiffness.Reluctivity(i).Full{j}{2,2}(K,J);              
%                    	
%                    	this.Stiffness.Reluctivity(i).Full{j}{1,1}(J,:) = [];
%                     this.Stiffness.Reluctivity(i).Full{j}{1,1}(:,J) = [];
%                     
%                   	this.Stiffness.Reluctivity(i).Full{j}{1,2}(J,:) = [];
%                     this.Stiffness.Reluctivity(i).Full{j}{1,2}(:,J) = [];
%                     
%                    	this.Stiffness.Reluctivity(i).Full{j}{2,1}(J,:) = [];
%                     this.Stiffness.Reluctivity(i).Full{j}{2,1}(:,J) = [];
%                     
%                   	this.Stiffness.Reluctivity(i).Full{j}{2,2}(J,:) = [];
%                     this.Stiffness.Reluctivity(i).Full{j}{2,2}(:,J) = [];
                                    
                    %% Mass
%                     warning('MotorProto:Verbose', 'Ensure proper behavior for conducting regions occuring on boundaries');
                   	
%                     this.Mass.Conductivity(i).Full{j}{1,1}(:) = 0;
%                     this.Mass.Conductivity(i).Full{j}{2,2}(:) = 0;
%                     this.Mass.Conductivity(i).Full{j}{1,2} ...
%                                 = -this.Mass.Conductivity(i).Full{j}{1,2};
%                     
%                    	this.Mass.Conductivity(i).Full{j}{1,1}(J,:) = [];
%                     this.Mass.Conductivity(i).Full{j}{1,1}(:,J) = [];
%                     
%                   	this.Mass.Conductivity(i).Full{j}{1,2}(J,:) = [];
%                     this.Mass.Conductivity(i).Full{j}{1,2}(:,J) = [];
%                     
%                    	this.Mass.Conductivity(i).Full{j}{2,1}(J,:) = [];
%                     this.Mass.Conductivity(i).Full{j}{2,1}(:,J) = [];
% 
%                   	this.Mass.Conductivity(i).Full{j}{2,2}(J,:) = [];
%                     this.Mass.Conductivity(i).Full{j}{2,2}(:,J) = []; 
                    
                    %% Curl dBydXz
%                     this.Jacobian.Magnetization(i).dBydXz{j}{1,1}(:,I)...
%                     	=            this.Jacobian.Magnetization(i).dBydXz{j}{1,1}(:,I) ...
%                           + C(1,1) * this.Jacobian.Magnetization(i).dBydXz{j}{1,1}(:,J);
%                       
%                    	this.Jacobian.Magnetization(i).dBydXz{j}{1,2}(:) = 0;
%                     this.Jacobian.Magnetization(i).dBydXz{j}{1,2}(:,I)...
%                     	= C(2,1) * this.Jacobian.Magnetization(i).dBydXz{j}{1,1}(:,J);
%                       
%                   	this.Jacobian.Magnetization(i).dBydXz{j}{2,2}(:,I)...
%                     	=            this.Jacobian.Magnetization(i).dBydXz{j}{2,2}(:,I) ...
%                           + C(2,2) * this.Jacobian.Magnetization(i).dBydXz{j}{2,2}(:,J);
%                       
%                    	this.Jacobian.Magnetization(i).dBydXz{j}{2,1}(:) = 0;
%                     this.Jacobian.Magnetization(i).dBydXz{j}{2,1}(:,I)...
%                     	= C(1,2) * this.Jacobian.Magnetization(i).dBydXz{j}{2,2}(:,J);
                    
                    %% Curl dBxdXz
%                     this.Jacobian.Magnetization(i).dBxdXz{j}{1,1}(:,I)...
%                     	=            this.Jacobian.Magnetization(i).dBxdXz{j}{1,1}(:,I) ...
%                           + C(1,1) * this.Jacobian.Magnetization(i).dBxdXz{j}{1,1}(:,J);
%                       
%                    	this.Jacobian.Magnetization(i).dBxdXz{j}{1,2}(:) = 0;
%                     this.Jacobian.Magnetization(i).dBxdXz{j}{1,2}(:,I)...
%                     	= C(2,1) * this.Jacobian.Magnetization(i).dBxdXz{j}{1,1}(:,J);
%                       
%                   	this.Jacobian.Magnetization(i).dBxdXz{j}{2,2}(:,I)...
%                     	=            this.Jacobian.Magnetization(i).dBxdXz{j}{2,2}(:,I) ...
%                           + C(2,2) * this.Jacobian.Magnetization(i).dBxdXz{j}{2,2}(:,J);
%                       
%                    	this.Jacobian.Magnetization(i).dBxdXz{j}{2,1}(:) = 0;
%                     this.Jacobian.Magnetization(i).dBxdXz{j}{2,1}(:,I)...
%                     	= C(1,2) * this.Jacobian.Magnetization(i).dBxdXz{j}{2,2}(:,J);
                    
                    %% Curl dIzdMy
%                     this.Jacobian.Magnetization(i).dIzdMy{j}{1,1}(I,:)...
%                     	=            this.Jacobian.Magnetization(i).dIzdMy{j}{1,1}(I,:) ...
%                           + C(1,1) * this.Jacobian.Magnetization(i).dIzdMy{j}{1,1}(J,:);
%                       
%                    	this.Jacobian.Magnetization(i).dIzdMy{j}{1,2}(:) = 0;
%                     this.Jacobian.Magnetization(i).dIzdMy{j}{1,2}(I,:)...
%                     	= C(1,2) * this.Jacobian.Magnetization(i).dIzdMy{j}{1,1}(J,:);
%                       
%                   	this.Jacobian.Magnetization(i).dIzdMy{j}{2,2}(I,:)...
%                     	=            this.Jacobian.Magnetization(i).dIzdMy{j}{2,2}(I,:) ...
%                           + C(2,2) * this.Jacobian.Magnetization(i).dIzdMy{j}{2,2}(J,:);
%                       
%                    	this.Jacobian.Magnetization(i).dIzdMy{j}{2,1}(:) = 0;
%                     this.Jacobian.Magnetization(i).dIzdMy{j}{2,1}(I,:)...
%                     	= C(2,1) * this.Jacobian.Magnetization(i).dIzdMy{j}{2,2}(J,:);
                    
                    %% Curl dIzdMx
%                     this.Jacobian.Magnetization(i).dIzdMx{j}{1,1}(I,:)...
%                     	=            this.Jacobian.Magnetization(i).dIzdMx{j}{1,1}(I,:) ...
%                           + C(1,1) * this.Jacobian.Magnetization(i).dIzdMx{j}{1,1}(J,:);
%                       
%                    	this.Jacobian.Magnetization(i).dIzdMx{j}{1,2}(:) = 0;
%                     this.Jacobian.Magnetization(i).dIzdMx{j}{1,2}(I,:)...
%                     	= C(1,2) * this.Jacobian.Magnetization(i).dIzdMx{j}{1,1}(J,:);
%                       
%                   	this.Jacobian.Magnetization(i).dIzdMx{j}{2,2}(I,:)...
%                     	=            this.Jacobian.Magnetization(i).dIzdMx{j}{2,2}(I,:) ...
%                           + C(2,2) * this.Jacobian.Magnetization(i).dIzdMx{j}{2,2}(J,:);
%                       
%                    	this.Jacobian.Magnetization(i).dIzdMx{j}{2,1}(:) = 0;
%                     this.Jacobian.Magnetization(i).dIzdMx{j}{2,1}(I,:)...
%                     	= C(2,1) * this.Jacobian.Magnetization(i).dIzdMx{j}{2,2}(J,:);
                    
%                   	this.Jacobian.Magnetization(i).dBydXz{j}{1,1}(:,J) = [];
%                     this.Jacobian.Magnetization(i).dBydXz{j}{1,2}(:,J) = [];
%                    	this.Jacobian.Magnetization(i).dBydXz{j}{2,1}(:,J) = [];
%                     this.Jacobian.Magnetization(i).dBydXz{j}{2,2}(:,J) = [];
%                     
%                   	this.Jacobian.Magnetization(i).dBxdXz{j}{1,1}(:,J) = [];
%                     this.Jacobian.Magnetization(i).dBxdXz{j}{1,2}(:,J) = [];
%                    	this.Jacobian.Magnetization(i).dBxdXz{j}{2,1}(:,J) = [];
%                     this.Jacobian.Magnetization(i).dBxdXz{j}{2,2}(:,J) = [];
                    
%                   	this.Jacobian.Magnetization(i).dIzdMy{j}{1,1}(J,:) = [];
%                   	this.Jacobian.Magnetization(i).dIzdMy{j}{1,2}(J,:) = [];
%                   	this.Jacobian.Magnetization(i).dIzdMy{j}{2,1}(J,:) = [];
%                   	this.Jacobian.Magnetization(i).dIzdMy{j}{2,2}(J,:) = [];
%                     
%                   	this.Jacobian.Magnetization(i).dIzdMx{j}{1,1}(J,:) = [];
%                   	this.Jacobian.Magnetization(i).dIzdMx{j}{1,2}(J,:) = [];
%                     this.Jacobian.Magnetization(i).dIzdMx{j}{2,1}(J,:) = [];
%                   	this.Jacobian.Magnetization(i).dIzdMx{j}{2,2}(J,:) = [];
                    
                    %% Exogenous
%                    	this.Exogenous.Magnetic(i).Full{j}{1,1}(I,:) ...
%                     	=            this.Exogenous.Magnetic(i).Full{j}{1,1}(I,:) ...
%                           + C(1,1) * this.Exogenous.Magnetic(i).Full{j}{1,1}(J,:);
%                       
%                    	this.Exogenous.Magnetic(i).Full{j}{1,2}(:) = 0;
%                    	this.Exogenous.Magnetic(i).Full{j}{1,2}(I,:) ...
%                     	= C(1,2) * this.Exogenous.Magnetic(i).Full{j}{1,1}(J,:);
%                       
%                     this.Exogenous.Magnetic(i).Full{j}{2,1}(:) = 0;
%                    	this.Exogenous.Magnetic(i).Full{j}{2,1}(I,:) ...
%                     	= C(2,1) * this.Exogenous.Magnetic(i).Full{j}{2,2}(J,:);
%                       
%                    	this.Exogenous.Magnetic(i).Full{j}{2,2}(I,:) ...
%                     	=            this.Exogenous.Magnetic(i).Full{j}{2,2}(I,:) ...
%                           + C(2,2) * this.Exogenous.Magnetic(i).Full{j}{2,2}(J,:);
%                       
%                  	this.Exogenous.Magnetic(i).Full{j}{1,1}(J,:) = [];
%                     this.Exogenous.Magnetic(i).Full{j}{1,2}(J,:) = [];  
%                     this.Exogenous.Magnetic(i).Full{j}{2,1}(J,:) = [];
%                     this.Exogenous.Magnetic(i).Full{j}{2,2}(J,:) = [];
                    
                   	%% Radial Boundary
%                     for k = 1:2
%                         %% Fourier Transform Matrix D
%                         this.Boundary.Radial(i).Radius(k).D{j}{1,1}(:,I)...
%                             =            this.Boundary.Radial(i).Radius(k).D{j}{1,1}(:,I)...
%                      	      + C(1,1) * this.Boundary.Radial(i).Radius(k).D{j}{1,1}(:,J);
%                         
%                         this.Boundary.Radial(i).Radius(k).D{j}{1,2}(:) = 0;
%                         this.Boundary.Radial(i).Radius(k).D{j}{1,2}(:,I)...
%                             = C(2,1) * this.Boundary.Radial(i).Radius(k).D{j}{1,1}(:,J);
%                        	
%                         this.Boundary.Radial(i).Radius(k).D{j}{2,1}(:) = 0;
%                         this.Boundary.Radial(i).Radius(k).D{j}{2,1}(:,I)...
%                             = C(1,2) * this.Boundary.Radial(i).Radius(k).D{j}{2,2}(:,J);
%                           
%                         this.Boundary.Radial(i).Radius(k).D{j}{2,2}(:,I)...
%                             =            this.Boundary.Radial(i).Radius(k).D{j}{2,2}(:,I)...
%                      	      + C(2,2) * this.Boundary.Radial(i).Radius(k).D{j}{2,2}(:,J);
%                         
%                         this.Boundary.Radial(i).Radius(k).D{j}{1,1}(:,J) = [];
%                         this.Boundary.Radial(i).Radius(k).D{j}{1,2}(:,J) = [];
%                         this.Boundary.Radial(i).Radius(k).D{j}{2,1}(:,J) = [];
%                         this.Boundary.Radial(i).Radius(k).D{j}{2,2}(:,J) = [];
%                         
%                         %% Boundary Integral Matrix S
%                         this.Boundary.Radial(i).Radius(k).S{j}{1,1}(I,:)...
%                         	=            this.Boundary.Radial(i).Radius(k).S{j}{1,1}(I,:)...
%                               + C(1,1) * this.Boundary.Radial(i).Radius(k).S{j}{1,1}(J,:);
%                         
%                         this.Boundary.Radial(i).Radius(k).S{j}{1,2}(:) = 0;
%                      	this.Boundary.Radial(i).Radius(k).S{j}{1,2}(I,:)...
%                         	= C(1,2) * this.Boundary.Radial(i).Radius(k).S{j}{1,1}(J,:);
%                         
%                         this.Boundary.Radial(i).Radius(k).S{j}{2,1}(:) = 0;
%                      	this.Boundary.Radial(i).Radius(k).S{j}{2,1}(I,:)...
%                         	= C(2,1) * this.Boundary.Radial(i).Radius(k).S{j}{2,2}(J,:);
%                         
%                         this.Boundary.Radial(i).Radius(k).S{j}{2,2}(I,:)...
%                         	=            this.Boundary.Radial(i).Radius(k).S{j}{2,2}(I,:)...
%                               + C(2,2) * this.Boundary.Radial(i).Radius(k).S{j}{2,2}(J,:);
% 
%                         this.Boundary.Radial(i).Radius(k).S{j}{1,1}(J,:) = [];
%                         this.Boundary.Radial(i).Radius(k).S{j}{1,2}(J,:) = [];
%                         this.Boundary.Radial(i).Radius(k).S{j}{2,1}(J,:) = [];
%                         this.Boundary.Radial(i).Radius(k).S{j}{2,2}(J,:) = [];
%                         
%                         %% Inverse Fourier Transform Matrix P
%                         this.Boundary.Radial(i).Radius(k).P{j}(J,:) = [];
%                     end
                
%                     %% Post Processing
%                     this.PostProcessing(i).Reduced2Full{j}{1,1}(J,:) ...
%                         = C(1,1) * this.PostProcessing(i).Reduced2Full{j}{1,1}(I,:);
%                     
%                     this.PostProcessing(i).Reduced2Full{j}{1,2}(:) = 0;
%                     this.PostProcessing(i).Reduced2Full{j}{1,2}(J,:) ...
%                         = C(2,1) * this.PostProcessing(i).Reduced2Full{j}{1,1}(I,:);
%                     
%                     this.PostProcessing(i).Reduced2Full{j}{2,1}(:) = 0;
%                     this.PostProcessing(i).Reduced2Full{j}{2,1}(J,:) ...
%                         = C(1,2) * this.PostProcessing(i).Reduced2Full{j}{2,2}(I,:);
%                     
%                     this.PostProcessing(i).Reduced2Full{j}{2,2}(J,:) ...
%                         = C(2,2) * this.PostProcessing(i).Reduced2Full{j}{2,2}(I,:);
%                     
%                     this.PostProcessing(i).Reduced2Full{j}{1,1}(:,J) = [];
%                     this.PostProcessing(i).Reduced2Full{j}{1,2}(:,J) = [];
%                     this.PostProcessing(i).Reduced2Full{j}{2,1}(:,J) = [];
%                     this.PostProcessing(i).Reduced2Full{j}{2,2}(:,J) = [];
%                 end
                
%                 this.PostProcessing(i).Full2Reduced(J,:) = [];
                
            	%% Indices
                warning('MotorProto:Verbose', 'generalize remove indices / nodes');
                index.Local(i).A(J)  = [];
                index.Global(i).A(J) = [];
                d = cumsum(diff(index.Local(i).A) - 1);
                d = [0,d] + min(index.Local(i).A) - 1;
                
                index.Local(i).A  = index.Local(i).A  - d;
                index.Global(i).A = index.Global(i).A - d;
                
                if i ~= 1
                    index.Global(i).A = index.Global(i).A ...
                      - (min(index.Global(i).A) - max(index.Global(i-1).A) - 1);
                end
                
                index.Local(i).X(J)  = [];
                index.Global(i).X(J) = [];
                d = cumsum(diff(index.Local(i).X) - 1);
                d = [0,d] + min(index.Local(i).X) - 1;
                
                index.Local(i).X  = index.Local(i).X  - d;
                index.Global(i).X = index.Global(i).X - d;
                
                if i ~= 1
                    index.Global(i).X = index.Global(i).X ...
                      - (min(index.Global(i).X) - max(index.Global(i-1).X) - 1);
                end
                
                index.Local(i).Unknowns  = index.Local(i).Unknowns  - numel(J);
                
                for j = i:nMesh
                    index.Global(j).Unknowns ...
                                          = index.Global(j).Unknowns - numel(J);
                end

                index.Local(i).Regions  = index.Local(i).Regions  - numel(J);
                index.Global(i).Regions = index.Global(i).Regions - numel(J);
            end
            this.Index = index;
        end
        
        function this = applyRadialBoundaryConditions(this)
            for i = 1:numel(this.Mesh_)
                for j = 1:numel(this.Stiffness.Reluctivity(i).Full)
%                     for k = 1:numel(this.Stiffness.Reluctivity(i).Full{j})
%                         this.Stiffness.Reluctivity(i).Full{j}{k} ...
%                             = this.Stiffness.Reluctivity(i).Full{j}{k} ...
%                                 + real(this.Boundary.Radial.S(i).Inner{j}{k} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{k}) ...
%                                 + real(this.Boundary.Radial.S(i).Outer{j}{k} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{k});
%                     end
                   this.Stiffness.Reluctivity(i).Full{j}{1,1} = this.Stiffness.Reluctivity(i).Full{j}{1,1} ...
                        + real(   this.Boundary.Radial.S(i).Inner{j}{1,1} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{1,1} ...
                                + this.Boundary.Radial.S(i).Inner{j}{1,2} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{2,1}) ...
                        + real(   this.Boundary.Radial.S(i).Outer{j}{1,1} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{1,1} ...
                                + this.Boundary.Radial.S(i).Outer{j}{1,2} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{2,1});
                            
                    this.Stiffness.Reluctivity(i).Full{j}{2,1} = this.Stiffness.Reluctivity(i).Full{j}{2,1} ...
                        + real(   this.Boundary.Radial.S(i).Inner{j}{2,1} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{1,1} ...
                                + this.Boundary.Radial.S(i).Inner{j}{2,2} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{2,1}) ...
                        + real(   this.Boundary.Radial.S(i).Outer{j}{2,1} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{1,1} ...
                                + this.Boundary.Radial.S(i).Outer{j}{2,2} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{2,1});
                            
                    this.Stiffness.Reluctivity(i).Full{j}{1,2} = this.Stiffness.Reluctivity(i).Full{j}{1,2} ...
                        + real(   this.Boundary.Radial.S(i).Inner{j}{1,1} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{1,2} ...
                                + this.Boundary.Radial.S(i).Inner{j}{1,2} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{2,2}) ...
                        + real(   this.Boundary.Radial.S(i).Outer{j}{1,1} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{1,2} ...
                                + this.Boundary.Radial.S(i).Outer{j}{1,2} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{2,2});
                            
                    this.Stiffness.Reluctivity(i).Full{j}{2,2} = this.Stiffness.Reluctivity(i).Full{j}{2,2} ...
                        + real(   this.Boundary.Radial.S(i).Inner{j}{2,1} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{1,2} ...
                                + this.Boundary.Radial.S(i).Inner{j}{2,2} * this.Boundary.Radial.F(i).Inner * this.Boundary.Radial.D(i).Inner{j}{2,2}) ...
                        + real(   this.Boundary.Radial.S(i).Outer{j}{2,1} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{1,2} ...
                                + this.Boundary.Radial.S(i).Outer{j}{2,2} * this.Boundary.Radial.F(i).Outer * this.Boundary.Radial.D(i).Outer{j}{2,2});
                end
            end
        end
        
        %% Matrix/Vector Functions
        function linearMatrix = K(this, ~, h)
            assembly     = this.Assemblies_;
            nAssemblies  = numel(assembly);
            linearMatrix = cell(nAssemblies, nAssemblies);
            Nh           = cellfun('length',h);
            for i = 1:nAssemblies
                linearMatrix{i,i} = cell(1,Nh(i));
                for j = 1:Nh(i)
                    k = mod(h{i}(j),6) + 1;
                    if h{i}(j) ~= 0
                        linearMatrix{i,i}{j} = 2 * cell2mat(this.Stiffness.Reluctivity(i).Full{k});
                    else
                        linearMatrix{i,i}{j} = this.Stiffness.Reluctivity(i).Full{k}{1,1};
                    end
                end
                linearMatrix{i,i} = blkdiag(linearMatrix{i,i}{:});
                
                if i > 1
                 	linearMatrix{i,i-1} = cell(Nh(i),Nh(i-1));
                    
                    %% calculate rotation frequencies
                    R  = - this.Boundary.Radial.R(i-1).Outer + this.Boundary.Radial.R(i).Inner;
                    m  = -round(R ./ min(abs(R)));
                    
                    hplus  = bsxfun(@plus,h{i}.',h{i-1});
                    hminus = bsxfun(@minus,h{i}.',h{i-1});
                    
                    Rcc = sparse([],[],[],length(m),length(m),4);
                    Rcc = repmat({Rcc},size(hplus));
                    Rcs = Rcc;
                    Rsc = Rcc;
                    Rss = Rcc;
                    
                    for j = 1:Nh(i)
                        for k = 1:Nh(i-1)
                            if h{i}(j) ~= 0 && h{i-1}(k) ~= 0
                                I = (m == hplus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                Rsc{j,k}(I,I) = 1i;
                                Rss{j,k}(I,I) = -1;

                                I = (m == -hplus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                Rsc{j,k}(I,I) = -1i;
                                Rss{j,k}(I,I) = -1;

                                I = (m == hminus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                Rsc{j,k}(I,I) = 1i;
                                Rss{j,k}(I,I) = 1;

                                I = (m == -hminus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                Rsc{j,k}(I,I) = -1i;
                                Rss{j,k}(I,I) = 1;
 
                            elseif h{i}(j) == 0 && h{i-1}(k) ~= 0
                                I = (m == h{i-1}(k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                
                                I = (m == -h{i-1}(k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                
                            elseif h{i}(j) ~= 0 && h{i-1}(k) == 0
                                I = (m == h{i}(j));
                                Rcc{j,k}(I,I) = 1;
                                Rsc{j,k}(I,I) = 1i;
                               
                                I = (m == -h{i}(j));
                                Rcc{j,k}(I,I) = 1;
                                Rsc{j,k}(I,I) = -1i;
                                
                            else
                                I = (m == 0);
                                Rcc{j,k}(I,I) = 1;
                            end
                        end
                    end
                    
                    S = this.Boundary.Radial.S(i).Inner;
                    D = this.Boundary.Radial.D(i-1).Outer;
                    G = this.Boundary.Radial.G(i-1).Outer;
                    
                    for j = 1:Nh(i)
                        jj = mod(h{i}(j),6) + 1;
                        for k = 1:Nh(i-1)
                            kk = mod(h{i-1}(k),6) + 1;
                            if h{i}(j) ~= 0 && h{i-1}(k) ~=0
                                linearMatrix{i,i-1}{j,k} ...
                                    = cell2mat(S{jj})  ...
                                        * [ Rcc{j,k}*G  Rcs{j,k}*G;
                                         	Rsc{j,k}*G  Rss{j,k}*G ] ...
                                        * cell2mat(D{kk});
                            elseif h{i}(j) == 0 && h{i-1}(k) ~= 0
                            	linearMatrix{i,i-1}{j,k} ...
                                    = cell2mat(S{jj}(1,1)) ...
                                        * [Rcc{j,k}*G,Rcs{j,k}*G;] ...
                                        * cell2mat(D{kk});
                            elseif h{i}(j) ~= 0 && h{i-1}(k) == 0
                                linearMatrix{i,i-1}{j,k} ...
                                    = cell2mat(S{jj}) ...
                                        * [ Rcc{j,k}*G;
                                         	Rsc{j,k}*G] ...
                                        * cell2mat(D{kk}(1,1));
                            else
                                linearMatrix{i,i-1}{j,k} ...
                                    = cell2mat(S{jj}(1,1))  ...
                                        * Rcc{j,k}*G ...
                                        * cell2mat(D{kk}(1,1));
                            end
                        end
                    end
                    linearMatrix{i,i-1} = real(cell2mat(linearMatrix{i,i-1}));
                end
                
                if i < nAssemblies
                    linearMatrix{i,i+1} = cell(Nh(i),Nh(i+1));
                    
                    %% calculate rotation frequencies
                    R  = - this.Boundary.Radial.R(i+1).Inner + this.Boundary.Radial.R(i).Outer;
                    m  = -round(R / min(abs(R)));
                    
                    hplus  = bsxfun(@plus,h{i}.',h{i+1});
                    hminus = bsxfun(@minus,h{i}.',h{i+1});
                    
                    Rcc = sparse([],[],[],length(m),length(m),4);
                    Rcc = repmat({Rcc},size(hplus));
                    Rcs = Rcc;
                    Rsc = Rcc;
                    Rss = Rcc;
                    
                    for j = 1:Nh(i)
                        for k = 1:Nh(i+1)
                            if h{i}(j) ~= 0 && h{i+1}(k) ~= 0
                                I = (m == hplus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                Rsc{j,k}(I,I) = 1i;
                                Rss{j,k}(I,I) = -1;

                                I = (m == -hplus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                Rsc{j,k}(I,I) = -1i;
                                Rss{j,k}(I,I) = -1;

                                I = (m == hminus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                Rsc{j,k}(I,I) = 1i;
                                Rss{j,k}(I,I) = 1;

                                I = (m == -hminus(j,k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                Rsc{j,k}(I,I) = -1i;
                                Rss{j,k}(I,I) = 1;

                            elseif h{i}(j) == 0 && h{i+1}(k) ~= 0
                                I = (m == h{i+1}(k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = 1i;
                                
                                I = (m == -h{i+1}(k));
                                Rcc{j,k}(I,I) = 1;
                                Rcs{j,k}(I,I) = -1i;
                                
                            elseif h{i}(j) ~= 0 && h{i+1}(k) == 0
                                I = (m == h{i}(j));
                                Rcc{j,k}(I,I) = 1;
                                Rsc{j,k}(I,I) = 1i;
                               
                                I = (m == -h{i}(j));
                                Rcc{j,k}(I,I) = 1;
                                Rsc{j,k}(I,I) = -1i;
                                
                            else
                                I = (m == 0);
                                Rcc{j,k}(I,I) = 1;
                            end
                        end
                    end
                    
                    S = this.Boundary.Radial.S(i).Outer;
                    D = this.Boundary.Radial.D(i+1).Inner;
                    G = this.Boundary.Radial.G(i+1).Inner;
                    
                    for j = 1:Nh(i)
                        jj = mod(h{i}(j),6) + 1;
                        for k = 1:Nh(i+1)
                            kk = mod(h{i+1}(k),6) + 1;
                            if h{i}(j) ~= 0 && h{i+1}(k) ~=0
                                linearMatrix{i,i+1}{j,k} ...
                                    = cell2mat(S{jj})  ...
                                        * [ Rcc{j,k}*G   Rcs{j,k}*G;
                                         	Rsc{j,k}*G   Rss{j,k}*G ] ...
                                        * cell2mat(D{kk});
                            elseif h{i}(j) == 0 && h{i+1}(k) ~=0
                            	linearMatrix{i,i+1}{j,k} ...
                                    = cell2mat(S{jj}(1,1)) ...
                                        * [Rcc{j,k}*G,Rcs{j,k}*G;] ...
                                        * cell2mat(D{kk});
                            elseif h{i}(j) ~= 0 && h{i+1}(k) ==0
                                linearMatrix{i,i+1}{j,k} ...
                                    = cell2mat(S{jj}) ...
                                        * [ Rcc{j,k}*G;
                                            Rsc{j,k}*G] ...
                                        * cell2mat(D{kk}(1,1));
                            else
                                linearMatrix{i,i+1}{j,k} ...
                                    = cell2mat(S{jj}(1,1))  ...
                                        * Rcc{j,k} * G ...
                                        * cell2mat(D{kk}(1,1));
                            end
                        end
                    end
                    linearMatrix{i,i+1} = real(cell2mat(linearMatrix{i,i+1}));
                end
            end
            linearMatrix = cell2mat(linearMatrix);
        end
        
        function conductivityMatrix                     = C(this, t, h)
            warning('MotorProto:Verbose','Ensure consistent electrical frequency definitions');
            nAssemblies        = numel(h);
            conductivityMatrix = cell(1,nAssemblies);
            conductivity       = this.Mass.Conductivity;
            T                  = t(2) + t(end);
            w                  = 2 * pi / T;
            
            for i = 1:nAssemblies
                Nh = numel(h{i});
                conductivityMatrix{i} = cell(1,Nh);
                
                k = mod(h{i}(1),6) + 1;
                if h{i}(1) == 0
                    conductivityMatrix{i}{1} ...
                        = 0 * cell2mat(conductivity(i).Full{k}(1,1));
                else
                    conductivityMatrix{i}{1} ...
                        = 2 * w * h{i}(1) * cell2mat(conductivity(i).Full{k});
                end
                
                for j = 2:Nh
                    k = mod(h{i}(1),6) + 1;
                    conductivityMatrix{i}{j} ...
                        = 2 * w * h{i}(j) * cell2mat(conductivity(i).Full{k});
                end
                
                conductivityMatrix{i} = blkdiag(conductivityMatrix{i}{:});
            end
            conductivityMatrix = blkdiag(conductivityMatrix{:});
        end
        
        function [nonlinearJacobian, nonlinearFunction] = G(this, x, t, h)
            mesh        = this.Mesh_;
            assembly    = [mesh.Assembly];
            nAssemblies = numel(assembly);
            curlE2N     = this.Jacobian.MagnetizationCurrent;
            curlN2E     = this.Jacobian.FluxDensity;
            curl        = this.Del.Curl;
            
            nonlinearFunction = cell(nAssemblies,1);
            nonlinearJacobian = cell(1,nAssemblies);

            %% Oversample by a factor of two, converting to time domain
            t      = linspace(0, t(2)+t(end), 10*numel(t)+1);
            t(end) = [];
            x      = doPostProcessing(this, x, t, h);
            for i = 1:nAssemblies
                %% Find zero harmonic
                Nh = numel(h{i});
                
                %% Calculate Flux Density
                A  = cell2mat(x(i,1:(end-1)));
                Bx = curl(i).Zn2Xe * A;
                By = curl(i).Zn2Ye * A;
                
                %% Preallocate
                NZ     = size(Bx);
                Mx     = zeros(NZ);
                My     = zeros(NZ);
                dMxdBx = zeros(NZ);
                dMydBy = zeros(NZ);
                dMydBx = zeros(NZ);
                dMxdBy = zeros(NZ);     
                
                %% Calculate magnetization in the time domain;
                materials  = [assembly(i).Regions.Material];
                nMaterials = numel(materials);
              	for j = 1:nMaterials
                    if ~materials(j).Linear
                        J = (mesh(i).ElementRegions == j);
                        [Mx(J,:), My(J,:), dMxdBx(J,:), dMydBy(J,:), dMydBx(J,:), dMxdBy(J,:)] = materials(j).vectorM(Bx(J,:), By(J,:));
%                         Mx(J,:) = Bx(J,:)*(1-1/1000)/mu_o;
%                         My(J,:) = By(J,:)*(1-1/1000)/mu_o;
%                         dMxdBx(J,:) = (1-1/1000)/mu_o;
%                         dMydBy(J,:) = (1-1/1000)/mu_o;
%                         dMydBx(J,:) = 0;
%                         dMxdBy(J,:) = 0;
                    end
                end
                
                %% Transform the back to the frequency domain
                [cosMx,sinMx] = dscfft(Mx,[],2);
                [cosMy,sinMy] = dscfft(My,[],2);

               	[cosdMxdBx,sindMxdBx] = dscfft(dMxdBx,[],2);
                [cosdMydBy,sindMydBy] = dscfft(dMydBy,[],2);
            	[cosdMydBx,sindMydBx] = dscfft(dMydBx,[],2);
                [cosdMxdBy,sindMxdBy] = dscfft(dMxdBy,[],2);
                
                [nRows,nCols] = size(cosdMxdBx);
                
                %% Extract the necessary magnetization components
              	j         = h{i} + 1;
                j(j == 1) = [];
                Mx        = [2*cosMx(:,j); 2*sinMx(:,j)];
                My        = [2*cosMy(:,j); 2*sinMy(:,j)];
                Mx        = reshape(Mx,[],1);
                My        = reshape(My,[],1);
                
                if h{i}(1) == 0
                    Mx = [cosMx(:,1);Mx];
                    My = [cosMy(:,1);My];
                end
                    
                %% Create block curl matrices in cell format
                j         = mod(h{i},6) + 1;
                curlZn2Xe = curlN2E(i).dBxdXz(j);
                curlZn2Ye = curlN2E(i).dBydXz(j);
                curlXe2Zn = curlE2N(i).dIzdMx(j);
                curlYe2Zn = curlE2N(i).dIzdMy(j);
                
                if h{i}(1) == 0;
                    curlZn2Xe{1} = curlZn2Xe{1}(1,1);
                    curlZn2Ye{1} = curlZn2Ye{1}(1,1);
                    curlXe2Zn{1} = curlXe2Zn{1}(1,1);
                    curlYe2Zn{1} = curlYe2Zn{1}(1,1);
                end
                
                %% Create differential magnetization matrices in block cell form
                fCell = @(x)(sparse(1:numel(x),1:numel(x),x));
                
                dMxdBx      = cell(Nh,Nh);
                [dMxdBx{:}] = deal(cell(2,2));               	
                cosdMxdBx   = mat2cell(cosdMxdBx, nRows, ones(1,nCols));
                sindMxdBx   = mat2cell(sindMxdBx, nRows, ones(1,nCols));
                cosdMxdBx   = cellfun(fCell,cosdMxdBx, 'UniformOutput', false);
                sindMxdBx   = cellfun(fCell,sindMxdBx, 'UniformOutput', false);
                
                dMydBx      = cell(Nh,Nh);
                [dMydBx{:}] = deal(cell(2,2));  
                cosdMydBx   = mat2cell(cosdMydBx ,nRows,ones(1,nCols));
                sindMydBx   = mat2cell(sindMydBx ,nRows,ones(1,nCols));
                cosdMydBx   = cellfun(fCell,cosdMydBx,'UniformOutput',false);
                sindMydBx   = cellfun(fCell,sindMydBx,'UniformOutput',false);
                
                dMxdBy     = cell(Nh,Nh);
                [dMxdBy{:}] = deal(cell(2,2));
                cosdMxdBy   = mat2cell(cosdMxdBy ,nRows,ones(1,nCols));
                sindMxdBy   = mat2cell(sindMxdBy ,nRows,ones(1,nCols)); 
                cosdMxdBy   = cellfun(fCell,cosdMxdBy,'UniformOutput',false);
                sindMxdBy   = cellfun(fCell,sindMxdBy,'UniformOutput',false);
                
                dMydBy      = cell(Nh,Nh);
                [dMydBy{:}] = deal(cell(2,2));
                cosdMydBy   = mat2cell(cosdMydBy ,nRows,ones(1,nCols));
                sindMydBy   = mat2cell(sindMydBy ,nRows,ones(1,nCols));
                cosdMydBy   = cellfun(fCell,cosdMydBy,'UniformOutput',false);
                sindMydBy   = cellfun(fCell,sindMydBy,'UniformOutput',false);
                
                for j = 1:Nh
                    for k = 1:Nh
                        if h{i}(j) ~= 0 && h{i}(k) ~= 0
                            m = [h{i}(j) + h{i}(k), h{i}(j) - h{i}(k)];
                            s = sign(m(2));
                            m = abs(m) + 1;

                            dMxdBx{j,k} = 2 * ...
                                [ cosdMxdBx{m(1)} +     cosdMxdBx{m(2)},...
                             	  sindMxdBx{m(1)} - s * sindMxdBx{m(2)};
                              	  sindMxdBx{m(1)} + s * sindMxdBx{m(2)},...
                               	 -cosdMxdBx{m(1)} +     cosdMxdBx{m(2)}];
                            
                            dMxdBy{j,k} = 2 * ...
                                [ cosdMxdBy{m(1)} +     cosdMxdBy{m(2)},...
                             	  sindMxdBy{m(1)} - s * sindMxdBy{m(2)};
                              	  sindMxdBy{m(1)} + s * sindMxdBy{m(2)},...
                               	 -cosdMxdBy{m(1)} +     cosdMxdBy{m(2)}];
                            
                            dMydBx{j,k} = 2 * ...
                                [ cosdMydBx{m(1)} +     cosdMydBx{m(2)},...
                             	  sindMydBx{m(1)} - s * sindMydBx{m(2)};
                              	  sindMydBx{m(1)} + s * sindMydBx{m(2)},...
                               	 -cosdMydBx{m(1)} +     cosdMydBx{m(2)}];
                            
                            dMydBy{j,k} = 2 * ...
                                [ cosdMydBy{m(1)} +     cosdMydBy{m(2)},...
                             	  sindMydBy{m(1)} - s * sindMydBy{m(2)};
                              	  sindMydBy{m(1)} + s * sindMydBy{m(2)},...
                               	 -cosdMydBy{m(1)} +     cosdMydBy{m(2)}];
                                                
                        elseif h{i}(j) == 0 && h{i}(k) ~= 0
                            m = h{i}(k) + 1;
                            dMxdBx{j,k} = 2 * [cosdMxdBx{m}, sindMxdBx{m}];
                            dMxdBy{j,k} = 2 * [cosdMxdBy{m}, sindMxdBy{m}];
                            dMydBx{j,k} = 2 * [cosdMydBx{m}, sindMydBx{m}];
                            dMydBy{j,k} = 2 * [cosdMydBy{m}, sindMydBy{m}];
                            
                        elseif h{i}(j) ~= 0 && h{i}(k) == 0
                            m = h{i}(j) + 1;
                            dMxdBx{j,k} = 2 * [cosdMxdBx{m}; sindMxdBx{m}];
                            dMxdBy{j,k} = 2 * [cosdMxdBy{m}; sindMxdBy{m}];
                            dMydBx{j,k} = 2 * [cosdMydBx{m}; sindMydBx{m}];
                            dMydBy{j,k} = 2 * [cosdMydBy{m}; sindMydBy{m}];
                            
                        else
                            dMxdBx{j,k} = cosdMxdBx{1};
                            dMxdBy{j,k} = cosdMxdBy{1};
                            dMydBx{j,k} = cosdMydBx{1};
                            dMydBy{j,k} = cosdMydBy{1};
                        end
                    end
                end
                
                %% Change from cell to matrix format
                fCell     = @(x)(cell2mat(x));
                curlZn2Xe = cellfun(fCell,curlZn2Xe, 'UniformOutput', false);
                curlZn2Ye = cellfun(fCell,curlZn2Ye, 'UniformOutput', false);
                curlXe2Zn = cellfun(fCell,curlXe2Zn, 'UniformOutput', false);
                curlYe2Zn = cellfun(fCell,curlYe2Zn, 'UniformOutput', false);
                
                curlZn2Xe = blkdiag(curlZn2Xe{:});
                curlZn2Ye = blkdiag(curlZn2Ye{:});
                curlXe2Zn = blkdiag(curlXe2Zn{:});
                curlYe2Zn = blkdiag(curlYe2Zn{:});
                
                dMxdBx    = cell2mat(dMxdBx);
                dMxdBy    = cell2mat(dMxdBy);
                dMydBx    = cell2mat(dMydBx);
                dMydBy    = cell2mat(dMydBy);
                
                %% Perform matrix multiplication
                nonlinearFunction{i} = curlXe2Zn * Mx + curlYe2Zn * My; 
                nonlinearJacobian{i} =   (curlXe2Zn * dMxdBx + curlYe2Zn * dMydBx) * curlZn2Xe ...
                                       + (curlXe2Zn * dMxdBy + curlYe2Zn * dMydBy) * curlZn2Ye;
            end
            
            nonlinearJacobian = blkdiag(nonlinearJacobian{:});
            nonlinearFunction = cell2mat(nonlinearFunction);
        end
        
     	function exogenousFunction                      = f(this, t, h)
            mesh     = this.Mesh_;
            assembly = [mesh.Assembly];
            
            nAssemblies       = numel(assembly);
            exogenousFunction = cell(nAssemblies,1);
            curl              = this.Jacobian.MagnetizationCurrent;
            
            for i = 1:nAssemblies
                Nh = numel(h{i});
                E  = cell(1,Nh);
                for j = 1:Nh
                 	k = mod(h{i}(j),6) + 1;
                    if h{i}(j) ~= 0
                        E{j} = cell2mat(this.Exogenous.Magnetic(i).Full{k});
                    else
                        E{j} = this.Exogenous.Magnetic(i).Full{k}{1,1};
                    end
                end
                E = blkdiag(E{:});
                
               	sources = assembly(i).Sources;
                if numel(sources) > 0
                    f = sources.f(t, h{i});
                    f = dscfft(f,[],2);
                    
                    j = h{i}(h{i} ~= 0);
                    j = [2*j; 2*j+1];
                    j = reshape(j,1,[]);
                    if any(h{i} == 0)
                        j = [1 j];
                    end
                    
                    f = f(:,j);
                    f = reshape(f,[],1);
                    f = 2 * E * f;
                else
                    f = sparse(max(size(E)),1);
                end
                
                isZero = (h{i} == 0);
                if any(isZero)
                    materials = [assembly(i).Regions.Material];
                    
                    [~,Mx,My] = materials.vectorMr;
                    Mx        = sparse(Mx(mesh(i).ElementRegions));
                    My        = sparse(My(mesh(i).ElementRegions));
                    
                    Im        = - curl(i).dIzdMx{1}{1,1} * Mx - curl(i).dIzdMy{1}{1,1} * My;
                    
                    fRows     = numel(f) / (2*Nh-1) * ones(1,(2*Nh-1));
                    f         = mat2cell(f,fRows,1);
                    f{1}      = f{isZero} + Im;
                    f         = cell2mat(f);
                end
                
                exogenousFunction{i,1} = f;
            end
            exogenousFunction = cell2mat(exogenousFunction);
        end
        
       	%% Postprocessing Functions
        function [y, y_t] = doPostProcessing(this, x, t, h)
            nTimes      = numel(t);
            ppMatrices  = this.PostProcessing;
            nAssemblies = numel(ppMatrices);
            y           = cell(nAssemblies, nTimes);
            y_t         = cell(nAssemblies, nTimes);
            r2f         = cell(1, nAssemblies);
            nRows       = zeros(1, nAssemblies);
            
            for i = 1:nAssemblies
                Nh = numel(h{i});
                r2f{i} = cell(1,Nh);
                k      = mod(h{i},6) + 1;
                r2f{i} = ppMatrices(i).Reduced2Full(k);
                
                if h{i}(1) == 0
                    r2f{i}{1} = r2f{i}{1}(1,1);
                end
                
                r2f{i} = cellfun(@(x)(cell2mat(x)),r2f{i},'UniformOutput',false);
                r2f{i} = blkdiag(r2f{i}{:});
                
                [nRows(i),~] = size(r2f{i});
            end
            
            r2f       = blkdiag(r2f{:});
            [~,nCols] = size(r2f);
            r2f       = mat2cell(r2f, nRows, nCols);
            
            for i = 1:nAssemblies
                if h{i}(1) == 0
                    Nh        = 2 * numel(h{i}) - 1;
                    I         = 2*h{i}(2:end);
                    I         = [I;I+1];
                    I         = [1, reshape(I,1,[])];
                else
                    Nh        = 2 * numel(h{i});
                    I         = 2 * h{i};
                    I         = [I;I+1];
                    I         = reshape(I,1,[]);
                end
                
                xh        = r2f{i} * x;
                xh        = reshape(xh, [], Nh);
                nRows(i)  = nRows(i) / Nh;
                
                yh        = zeros(nRows(i), nTimes);
                yh(:,I)   = xh;
                
                if mod(nTimes,2) == 0
                    Nh   = nTimes / 2 - 1;
                    NMax = nTimes - 1;
                else
                    Nh   = (nTimes - 1) / 2;
                    NMax = nTimes;
                end
                
                T                = t(2) + t(end);
                Omega            = 2 * pi / T * diag(1:Nh);
                yh_t             = yh;
                yh_t(:,2:2:NMax) = - yh(:,3:2:NMax) * Omega;
                yh_t(:,3:2:NMax) =   yh(:,2:2:NMax) * Omega;
                yh_t(:,1)        = 0;
                
                yh        = idscfft(yh,   [], 2);
                yh_t      = idscfft(yh_t, [], 2);
                
                y(i,:)    = mat2cell(yh,   nRows(i), ones(1,nTimes));
                y_t(i,:)  = mat2cell(yh_t, nRows(i), ones(1,nTimes));
            end
            
            y(:,end+1)   = y(:,1);
            y_t(:,end+1) = y_t(:,1);
        end
    end
end