classdef Circuit < Component
    properties
        ConnectionMatrices = {}
        ConnectionPolarity = {}
    end
    
    properties (Abstract, SetAccess = private)
        TerminalNames
    end
    
    methods
        function this = Circuit(varargin)
            this = this@Component(varargin{:});
        end
    end
end