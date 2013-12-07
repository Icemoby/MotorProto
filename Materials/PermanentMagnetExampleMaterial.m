classdef PermanentMagnetExampleMaterial < LinearMaterial
    properties (SetAccess = protected)
        Description         = 'Permanent Magnet';
        Color               = 'm';
        
        RemanentFluxDensity = 1.2;
        Conductivity        = 1.82e6;
        Permeability        = mu_o;
    end
end