classdef NdFe35 < LinearMaterial
    properties (SetAccess = protected)
        Description         = 'NdFe35';
        Color               = 'm';
        
        RemanentFluxDensity = 1.23;
        Conductivity        = 1.82e6;
        Permeability        = mu_o * 1.09978;
    end
end