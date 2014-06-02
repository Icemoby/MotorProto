classdef FieldWoundTransformer < Circuit
    properties
        FieldSlots
        TransformerSlots
    end
    
    properties (SetAccess = private)
        TerminalNames = {'Field Winding', 'Transformer Winding'};
    end
    
    methods
        function this = FieldWoundTransformer(varargin)
            this = this@Circuit(varargin{:});
        end
        
        function [rows, cols, vals] = circuitMatrix(this)
            connectionMatrices = this.ConnectionMatrices;
            
            %% Pre-Allocation Phase
            nEntries = 0;
            for j = 1:2
                [~, nParallel] = size(connectionMatrices{j});
                nEntries       = nEntries + nParallel * 3;
            end

            rows = zeros(1, nEntries);
            cols = zeros(1, nEntries);
            vals = zeros(1, nEntries);
            
            %% Build Phase
            ind = 0;
            pos = 0;
            
            for j = 1:2
                [~, nParallel] = size(connectionMatrices{j});

                %FEA Voltage Equation
                rows(ind+1:ind+nParallel) = pos + (1:nParallel);
                cols(ind+1:ind+nParallel) = pos + nParallel + 1;
                vals(ind+1:ind+nParallel) = -1;

                ind = ind + nParallel;

                %Line Current Equation
                %Sum of currents in parallel paths in the finite element model equals the applied line currents
                rows(ind+1:ind+nParallel) = pos + nParallel + 1;
                cols(ind+1:ind+nParallel) = pos + (1:nParallel);
                vals(ind+1:ind+nParallel) = -1;

                ind = ind + nParallel;

                pos = pos + nParallel + 1;
            end
            
            %% Delete extra allocation space
            rows(ind+1:end) = [];
            cols(ind+1:end) = [];
            vals(ind+1:end) = [];
        end
        
        function [rows, cols, vals] = scalarPotential2CircuitMatrix(this)       
            connectionMatrices = this.ConnectionMatrices;
            connectionPolarity = this.ConnectionPolarity;

            %% Preallocation Phase
            nEntries = 0;
            for j = 1:2
                nEntries = nEntries + numel(connectionMatrices{j});
            end

            rows = zeros(1, nEntries);
            cols = zeros(1, nEntries);
            vals = zeros(1, nEntries);
            
            %% Build Phase
            ind = 0;
            pos = 0;   
            for j = 1:2
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
    end
    
    methods (Static)
        function componentOut = newComponent(varargin)
            componentOut = FieldWoundTransformer(varargin{:});
        end
    end
end