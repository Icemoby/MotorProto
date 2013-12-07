classdef Air < LinearMaterial
    properties (SetAccess = protected)
        Description         = 'Air';
        Color               = 'w';
        
        RemanentFluxDensity = 0;
        Permeability        = mu_o;
        Conductivity        = 0;
    end
end