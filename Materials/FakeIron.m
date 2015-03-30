classdef FakeIron < NonlinearMaterial
    properties (SetAccess = protected)
        %% Basic material information
        Description  = 'Arnon7';
        Color = 'b';

        %% BH Curve Data
        HData = [0, 10, 100, 1000];
        BData = [0, 10, 100, 1000] / (1000*mu_o);
      	
        BasisDegree = 2;
              
        %% Empirical Loss Data
        CoreLossData = [[0.015 0.035 0.600 0.950 3.000 0.080 0.400 1.500 3.000 5.000 7.000 0.30 0.90 4.00 8.00 1.50 3.50] * 1.7e4;
                        [200.0 200.0 200.0 200.0 200.0 400.0 400.0 400.0 400.0 400.0 1000  2000 2000 2000 2000 5000 5000];
                        [0.050 0.080 0.400 0.500 1.000 0.080 0.200 0.400 0.600 0.800 0.500 0.05 0.09 0.20 0.30 0.06 0.09]];

       
        RemanentFluxDensity = 0;
        Conductivity = 2.12e6;
        Density = 7650;
    end
end