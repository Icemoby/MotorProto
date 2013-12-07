classdef CurrentSource < Source
    properties (Constant)
        Type = SourceTypes.CurrentSource;
    end
    
    properties
    	ExternalResistance  = Source.setProperty(inf);
    end
    
    methods
        function this = CurrentSource(varargin)
            this                = this@Source(varargin{:});
            this.ConnectionType = ConnectionTypes.Delta;
        end
        
        function this = set.ExternalResistance(this, value)
            this.ExternalResistance = this.setProperty(value);
        end
        
    	function [rows, cols, vals] = circuitMatrix(this, connectionType)       
            connectionMatrices = this.ConnectionMatrices;
            nPhases            = this.Phases.Value;
            
            %% Pre-Allocation Phase
            nEntries           = nPhases * 6 + 1;
            for j = 1:nPhases
                [~, nParallel] = size(connectionMatrices{j});
                nEntries       = nEntries + nParallel * 3;
            end

            rows = zeros(1, nEntries);
            cols = zeros(1, nEntries);
            vals = zeros(1, nEntries);
            
            cmInd = zeros(1, nPhases);
            cmVal = zeros(1, nPhases);
            
            %% Build Phase
            ind = 0;
            pos = 0;
            for j = 1:nPhases
                [~, nParallel]            = size(connectionMatrices{j});
                
                %% FEA Voltage Equation
                rows(ind+1:ind+nParallel) = pos + (1:nParallel);
                cols(ind+1:ind+nParallel) = pos + nParallel + 1;
             	vals(ind+1:ind+nParallel) = -1;
                    
                ind = ind + nParallel;
                
                switch connectionType
                    case ConnectionTypes.Wye
                        %% Line Current Equation
                        %Sum of currents in parallel paths in the finite element model equals the applied line currents
                        rows(ind+1:ind+nParallel) = pos + nParallel + 1;
                        cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                        vals(ind+1:ind+nParallel) = -1;

                        ind = ind + nParallel;
                    case ConnectionTypes.Delta
                        %% Line Current Equation
                        %Sum of currents in parallel paths in the finite element model plus 1/3 FEA common mode current
                        %equals the applied phase current minus 1/3 source common mode current
                        rows(ind+1:ind+nParallel) = pos + nParallel + 1;
                        cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                        vals(ind+1:ind+nParallel) = -1;

                        ind = ind + nParallel;
                        
                        cmInd(j) = pos + nParallel + 1;
                        cmVal(j) = 1 / nPhases;
                    otherwise
                        error('MotorProto:VoltageSource','Unknown ConnectionTypes %s', char(connectionType))
                end
                pos = pos + nParallel + 1;
            end
            
            if connectionType == ConnectionTypes.Delta
                %% Add constraints/common mode values
                rows(ind+1:ind+nPhases) = cmInd;
                cols(ind+1:ind+nPhases) = pos + 1;
                vals(ind+1:ind+nPhases) = cmVal;

                ind = ind + nPhases;

                rows(ind+1:ind+nPhases) = pos + 1;
                cols(ind+1:ind+nPhases) = cmInd;
                vals(ind+1:ind+nPhases) = cmVal;

                ind = ind + nPhases;
            end
            
            rows(ind+1:end) = [];
            cols(ind+1:end) = [];
            vals(ind+1:end) = [];
        end
        
        function [rows, cols, vals] = scalarPotential2CircuitMatrix(this)       
            connectionMatrices = this.ConnectionMatrices;
            connectionPolarity = this.ConnectionPolarity;
            nPhases            = this.Phases.Value;

            %% Preallocation Phase
            nEntries = 0;
            for j = 1:nPhases
                nEntries = nEntries + numel(connectionMatrices{j});
            end

            rows = zeros(1, nEntries);
            cols = zeros(1, nEntries);
            vals = zeros(1, nEntries);
            
            %% Build Phase
            ind = 0;
            pos = 0;
            for j = 1:nPhases
                [nSeries, nParallel] = size(connectionMatrices{j});
                for k = 1:nParallel
                    rows(ind+1:ind+nSeries) = pos + k;
                    cols(ind+1:ind+nSeries) = connectionMatrices{j}(:,k);
                    vals(ind+1:ind+nSeries) = connectionPolarity{j}(:,k);
                    
                    ind = ind + nSeries;
                end
                pos = pos + nParallel + 1;
            end
        end
        
        function [rows, cols, vals] = source2CircuitMatrix(this, connectionType)
            connectionMatrices = this.ConnectionMatrices;
            nPhases            = this.Phases.Value;
            
            switch connectionType
                case ConnectionTypes.Delta
                    rows            = zeros(1,nPhases^2);
                    cols            = repmat(1:nPhases, 1, nPhases);
                    vals            = zeros(1,nPhases^2);
                    
                    [~, nParallel]  = size(connectionMatrices{1});
                    rows(1:nPhases) = nParallel + 1;
                    vals(1:nPhases) = [-(nPhases - 1) ones(1, nPhases - 1)] / nPhases;
                    for j = 2:nPhases
                        [~, nParallel] = size(connectionMatrices{j});
                        J              = ((j-1)*nPhases+1):(j*nPhases);
                        rows(J)        = rows(j-1) + nParallel + 2;
                        vals(J)        = [ones(1, j - 1) -(nPhases - 1) ones(1, nPhases - j)] / nPhases;
                    end
                case ConnectionTypes.Wye
                    rows          = zeros(1,2*nPhases);
                    cols          = [1:nPhases, circshift(1:nPhases, [1,1])];
                    vals          = [-ones(1, nPhases) ones(1, nPhases)];
                    
                    [~, nParallel] = size(connectionMatrices{1});
                    rows(1)        = nParallel + 1;
                    for j = 2:nPhases
                        [~, nParallel] = size(connectionMatrices{j});
                        rows(j)        = rows(j-1) + nParallel + 1;
                    end
                    rows(nPhases+1:end) = rows(1:nPhases);
                otherwise
                    error('MotorProto:VoltageSource','Unknown ConnectionTypes %s', char(connectionType))
            end
        end
    end
    
    methods (Static)
        function componentOut = newComponent(varargin)
            componentOut = CurrentSource(varargin{:});
        end
    end
end