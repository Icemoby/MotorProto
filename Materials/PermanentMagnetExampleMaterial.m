classdef PermanentMagnetExampleMaterial < LinearMaterial
    properties (SetAccess = protected)
        Description         = 'Permanent Magnet';
        Color               = 'm';
        
        RemanentFluxDensity = 1.2;
        Conductivity        = 6.67e5;
        Density             = 7500;
        Permeability        = mu_o;
    end
end