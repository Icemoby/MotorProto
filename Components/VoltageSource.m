classdef VoltageSource < Source
    properties (Constant)
        Type = SourceTypes.VoltageSource
    end
    
    properties
    	ExternalResistance  = Source.setProperty(0);
    end
    
    methods
        function this = VoltageSource(varargin)
            this                = this@Source(varargin{:});
            this.ConnectionType = ConnectionTypes.Wye;
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
                
                rows(ind+1:ind+nParallel) = pos + (1:nParallel);
                cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                vals(ind+1:ind+nParallel) = -0;
                
                ind = ind + nParallel;
                
                switch connectionType
                    case ConnectionTypes.Wye
                        %% Line Current Equation
                        %Sum of currents in parallel paths in the finite element model...
                        rows(ind+1:ind+nParallel) = pos + nParallel + 1;
                        cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                        vals(ind+1:ind+nParallel) = -1;

                        ind = ind + nParallel;
                        
                        %...equals the negative of the line currents
                        rows(ind+1)               = pos + nParallel + 1;
                        cols(ind+1)               = pos + nParallel + 2;
                        vals(ind+1)               = -1;

                        ind = ind + 1;
                        
                        %% Line Voltage Equation
                        %The finite element voltages...
                        rows(ind+1)               = pos + nParallel + 2;
                        cols(ind+1)               = pos + nParallel + 1;
                        vals(ind+1)               = -1;
                        
                        ind = ind + 1;
                        
                        %...plus the common mode voltage equals the line to neutral voltage
                        %Line currents must sum to zero
                        cmInd(j)                   = pos + nParallel + 2;
                        cmVal(j)                   = -1;
                    case ConnectionTypes.Delta
                        %% Line Current Equation
                        %The sum of currents in parallel paths in the finite element model...
                        rows(ind+1:ind+nParallel) = pos + nParallel + 1;
                        cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                        vals(ind+1:ind+nParallel) = -1;

                        ind = ind + nParallel;

                        %...equals the phase currents
                        rows(ind+1)               = pos + nParallel + 1;
                        cols(ind+1)               = pos + nParallel + 2;
                        vals(ind+1)               = 1;

                        ind = ind + 1;
                        
                        %% Line Voltage Equation
                        %The finite element voltages equal the line to line voltages
                        rows(ind+1)               = pos + nParallel + 2;
                        cols(ind+1)               = pos + nParallel + 1;
                        vals(ind+1)               = 1;
                        
                        ind = ind + 1;
                    otherwise
                        error('MotorProto:VoltageSource','Unknown ConnectionTypes %s', char(connectionType))
                end
                pos = pos + nParallel + 2;
            end
            
            if connectionType == ConnectionTypes.Wye
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
                pos = pos + nParallel + 2;
            end
        end
        
        function [rows, cols, vals] = source2CircuitMatrix(this, connectionType)
            connectionMatrices = this.ConnectionMatrices;
            nPhases            = this.Phases.Value;
            
            switch connectionType
                case ConnectionTypes.Wye
                    rows          = zeros(1,nPhases);
                    cols          = 1:nPhases;
                    vals          = ones(1, nPhases);
                    
                    [~,nParallel] = size(connectionMatrices{1});
                    rows(1)       = nParallel + 2;
                    for j = 2:nPhases
                        [~,nParallel] = size(connectionMatrices{j});
                        rows(j)       = rows(j-1) + nParallel + 2;
                    end
                case ConnectionTypes.Delta
                    rows          = zeros(1,2*nPhases);
                    cols          = [1:nPhases, circshift(1:nPhases, [1,1])];
                    vals          = [ones(1, nPhases) -ones(1, nPhases)];
                    
                    [~,nParallel] = size(connectionMatrices{1});
                    rows(1)       = nParallel + 2;
                    for j = 2:nPhases
                        [~,nParallel] = size(connectionMatrices{j});
                        rows(j)       = rows(j-1) + nParallel + 2;
                    end
                    rows(nPhases+1:end) = rows(1:nPhases);
                otherwise
                    error('MotorProto:VoltageSource','Unknown ConnectionTypes %s', char(connectionType))
            end
        end
    end
    
    methods (Static)
        function componentOut = newComponent(varargin)
            componentOut = VoltageSource(varargin{:});
        end
    end
end