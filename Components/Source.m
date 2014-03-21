classdef Source < Component
    properties
        ElectricalFrequency = Source.setProperty(60);      
        HarmonicNumbers     = Source.setProperty(1);
        HarmonicAmplitudes  = Source.setProperty(0);
        HarmonicPhases      = Source.setProperty(0);
        Phases              = Source.setProperty(3);
        
        ConnectionType      = ConnectionTypes.Wye;
        ConnectionMatrices  = {};
        ConnectionPolarity  = {};
    end
    
    properties (Abstract, Constant)
        Type
    end
    
    methods
        %% Constructor
        function this = Source(varargin)
            if nargin ~= 0
                if isa(varargin{1}, 'Source')
                    inputObject              = varargin{1};
                    this.Name                = inputObject.Name;
                    this.ElectricalFrequency = inputObject.ElectricalFrequency;
                    this.HarmonicAmplitudes  = inputObject.HarmonicAmplitudes;
                    this.HarmonicPhases      = inputObject.HarmonicPhases;
                    this.ExternalResistance  = inputObject.ExternalResistance;
                    this.Phases              = inputObject.Phases;
                    this.ConnectionType      = inputObject.ConnectionType;
                    this.ConnectionMatrices  = inputObject.ConnectionMatrices;
                    this.ConnectionPolarity  = inputObject.ConnectionPolarity;
                    for i = 2:2:(nargin-1)
                        this.(varargin{i}) = varargin{i+1};
                    end
                else
                    this.Name = varargin{1};
                    for i = 2:2:(nargin-1)
                        this.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        %% Setters
        function this = set.ElectricalFrequency(this, value)
            this.ElectricalFrequency = this.setProperty(value);
        end
        
        function this = set.HarmonicNumbers(this, value)
            this.HarmonicNumbers = this.setProperty(value);
        end
        
        function this = set.HarmonicAmplitudes(this, value)
            this.HarmonicAmplitudes = this.setProperty(value);
        end
        
        function this = set.HarmonicPhases(this, value)
            this.HarmonicPhases = this.setProperty(value);
        end
        
        function this = set.Phases(this, value)
            this.Phases = this.setProperty(value);
        end
        
        function this = set.ConnectionType(this, value)
            if ischar(value)
                this.ConnectionType = ConnectionTypes.(value);
            else
                this.ConnectionType = ConnectionTypes(value);
            end
        end
    end
    
    methods (Sealed)        
        function waveform = f(this, t, h)
            m        = this.Phases.Value;
            waveform = zeros(m, numel(t));
            
            if nargin == 2
                n = this.HarmonicNumbers.Value;
                J = true(size(n));
            else
                J = ismember(this.HarmonicNumbers.Value, h);
                n = this.HarmonicNumbers.Value(J);
            end
            
            w = 2 * pi * this.ElectricalFrequency.Value;
            A = this.HarmonicAmplitudes.Value(:,J);
            p = this.HarmonicPhases.Value(:,J);
            s = -(0:(m-1)) * 2 * pi / m;
            I = ones(1,numel(t));
            
            if numel(A) == 1
                for i = 1:m
                    waveform(i,:) = A * cos(w * n' * t + (p' - n' * s(i)) * I);
                end
            else
                for i = 1:m
                    waveform(i,:) = A(i,:) * cos(w * n' * t + (p(i,:)' - n' * s(i)) * I);
                end
            end
        end
    end
end