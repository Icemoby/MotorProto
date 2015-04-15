classdef NonlinearMaterial < MaterialProperty        
    %NonlinearMaterial.m A class representing materials with nonlinear magnetic properties
    %
    % NonlinearMaterial properties:
    %   HData                 - Magnetic field intensity data from the B-H curve 
    %   BData                 - Magnetic flux density data from the B-H curve 
    %   MData                 - Magnetization calculated from HData and BData
    %   CoreLossData          - LossDensity/Frequency/FluxDensity core-loss measurements
    %   CoreLossCoefficients  - Coefficients for the Generalized Steinmetz equation fit to the CoreLossData
    %
    % NonlinearMaterial methods:
    %   calculateCoreLossCoefficients - Calculates the CoreLossCoefficients based on the CoreLossData
    %
    %   NonlinearMaterial inherets properties and methods from
    %   MaterialProperty. See the help for MaterialProperty for more
    %   information.
    %
    % See also MaterialProperty

%{
properties
    %HData - 
    %
    % See also NonlinearMaterial
    HData;
    
    %BData - 
    %
    % See also NonlinearMaterial
    BData;
    
    %MData - 
    %
    % See also NonlinearMaterial
    MData;
    
    %CoreLossData - 
    %
    % See also NonlinearMaterial
    CoreLossData;
    
    %CoreLossCoefficients - 
    %
    % See also NonlinearMaterial
    CoreLossCoefficients;
%}

    properties (Abstract, SetAccess = protected)
        %% Magnetization Curve Data
        HData
        BData
        
        %% Empirical Loss Data
        CoreLossData
    end
    
    properties (SetAccess = protected)
        Linear = false;
        
        %% Nonlinear Data
        MData
        
        %% Spline Data
        Bs
        Chi0
        a
        
        %% Core Losses
        CoreLossCoefficients
    end
    
    methods 
        %% Constructor Methods
        function this = NonlinearMaterial(varargin)
            this = this@MaterialProperty(varargin{:});
            
            %% B-H Curve Preprocessing
            %   Ensure monotonicity of M-B curve which implies, but is not implied by,
            %   monotonicity of the B-H curve. M = B/mu_o - H.
            %
            %   Ensure that dMdB = 0 as B->infinity
            
            B = this.BData;
            if size(B,1) == 1
                B = B';
            end
            
            H = this.HData;
            if size(H,1) == 1
                H = H';
            end
            
            M = B / mu_o - H;
            this.MData = M;
            
            M(1) = [];
            B(1) = [];
            H(1) = [];
            
            K     = length(M);
            Bsv   = min(1.0,max(B));
            Chi0v = (M(2)-M(1)) / (B(2) - B(1));
            av    = 2;
            tol   = inf;
            
            while tol > sqrt(eps / 3)
                f    = 0;
                dm   = zeros(3,1);
                dmm  = zeros(3,3);

                g = zeros(3,1);
                h = zeros(3,3);

                for k = 1:K
                    x = B(k);
                    
                  	%ek = (Chi0v*(x/2 - (atan(av*(Bsv - x))*(Bsv - x))/pi) + (Bsv*Chi0v*atan(Bsv*av))/pi - M(k));
                    %f = f + 0.5*ek^2 / H(k)^2;

                    %ek = ((Chi0v*x)/2 + (Bsv*Chi0v*erf(Bsv*av))/2 - (Chi0v*erf(av*(Bsv - x))*(Bsv - x))/2 + Chi0v/(2*pi^(1/2)*av*exp(Bsv^2*av^2)) - Chi0v/(2*pi^(1/2)*av*exp(av^2*(Bsv - x)^2))) - (M(k));
                    ek = log((Chi0v*x)/2 + (Bsv*Chi0v*erf(Bsv*av))/2 - (Chi0v*erf(av*(Bsv - x))*(Bsv - x))/2 + Chi0v/(2*pi^(1/2)*av*exp(Bsv^2*av^2)) - Chi0v/(2*pi^(1/2)*av*exp(av^2*(Bsv - x)^2))) - log(M(k));
                    f = f + 0.5*ek^2 / log(H(k))^2;
                    
                    m = (Chi0v.*x)./2 + (Bsv.*Chi0v.*erf(Bsv.*av))./2 - (Chi0v.*erf(av.*(Bsv - x)).*(Bsv - x))./2 + Chi0v./(2.*pi.^(1./2).*av.*exp(Bsv.^2.*av.^2)) - Chi0v./(2.*pi.^(1./2).*av.*exp(av.^2.*(Bsv - x).^2));
                    
                    dm(1) = (Chi0v*erf(Bsv*av))/2 - (Chi0v*erf(av*(Bsv - x)))/2 - (Chi0v*av*(Bsv - x))/(pi^(1/2)*exp(av^2*(Bsv - x)^2)) + (Chi0v*av*(2*Bsv - 2*x))/(2*pi^(1/2)*exp(av^2*(Bsv - x)^2));
                    dm(2) = x/2 + (Bsv*erf(Bsv*av))/2 - (erf(av*(Bsv - x))*(Bsv - x))/2 - 1/(2*pi^(1/2)*av*exp(av^2*(Bsv - x)^2)) + 1/(2*pi^(1/2)*av*exp(Bsv^2*av^2));
                    dm(3) = Chi0v/(2*pi^(1/2)*av^2*exp(av^2*(Bsv - x)^2)) - Chi0v/(2*pi^(1/2)*av^2*exp(Bsv^2*av^2));
                    
                    dmm(1,1) = (Chi0v*av)/(pi^(1/2)*exp(Bsv^2*av^2)) - (Chi0v*av)/(pi^(1/2)*exp(av^2*(Bsv - x)^2)) - (Chi0v*av^3*(2*Bsv - 2*x)^2)/(2*pi^(1/2)*exp(av^2*(Bsv - x)^2)) + (Chi0v*av^3*(Bsv - x)*(2*Bsv - 2*x))/(pi^(1/2)*exp(av^2*(Bsv - x)^2));
                    dmm(1,2) = erf(Bsv*av)/2 - erf(av*(Bsv - x))/2 - (av*(Bsv - x))/(pi^(1/2)*exp(av^2*(Bsv - x)^2)) + (av*(2*Bsv - 2*x))/(2*pi^(1/2)*exp(av^2*(Bsv - x)^2));
                    dmm(1,3) = (Chi0v*(2*Bsv - 2*x))/(2*pi^(1/2)*exp(av^2*(Bsv - x)^2)) - (2*Chi0v*(Bsv - x))/(pi^(1/2)*exp(av^2*(Bsv - x)^2)) + (Bsv*Chi0v)/(pi^(1/2)*exp(Bsv^2*av^2)) + (2*Chi0v*av^2*(Bsv - x)^3)/(pi^(1/2)*exp(av^2*(Bsv - x)^2)) - (Chi0v*av^2*(Bsv - x)^2*(2*Bsv - 2*x))/(pi^(1/2)*exp(av^2*(Bsv - x)^2));
                    dmm(2,2) = 0;
                    dmm(2,3) = 1/(2*pi^(1/2)*av^2*exp(av^2*(Bsv - x)^2)) - 1/(2*pi^(1/2)*av^2*exp(Bsv^2*av^2));
                    dmm(3,3) = Chi0v/(pi^(1/2)*av^3*exp(Bsv^2*av^2)) - Chi0v/(pi^(1/2)*av^3*exp(av^2*(Bsv - x)^2)) - (Chi0v*(Bsv - x)^2)/(pi^(1/2)*av*exp(av^2*(Bsv - x)^2)) + (Bsv^2*Chi0v)/(pi^(1/2)*av*exp(Bsv^2*av^2));
                    
%                     dm(1) = (Chi0v*atan(Bsv*av))/pi - Chi0v*(atan(av*(Bsv - x))/pi + (av*(Bsv - x))/(pi*(av^2*(Bsv - x)^2 + 1))) + (Bsv*Chi0v*av)/(pi*(Bsv^2*av^2 + 1));
%                     dm(2) = x/2 - (atan(av*(Bsv - x))*(Bsv - x))/pi + (Bsv*atan(Bsv*av))/pi;
%                     dm(3) = (Bsv^2*Chi0v)/(pi*(Bsv^2*av^2 + 1)) - (Chi0v*(Bsv - x)^2)/(pi*(av^2*(Bsv - x)^2 + 1));
% 
%                     dmm(1,1) = -(2*Chi0v*av^3*x*(2*Bsv - x)*(2*Bsv^2*av^2 - 2*Bsv*av^2*x + av^2*x^2 + 2))/(pi*(Bsv^2*av^2 + 1)^2*(Bsv^2*av^2 - 2*Bsv*av^2*x + av^2*x^2 + 1)^2);
%                     dmm(1,2) = atan(Bsv*av)/pi - atan(av*(Bsv - x))/pi + (Bsv*av)/(pi*(Bsv^2*av^2 + 1)) - (av*(Bsv - x))/(pi*(av^2*(Bsv - x)^2 + 1));
%                     dmm(1,3) = (2*Chi0v*x*(- 3*Bsv^4 + 6*Bsv^3*x - 4*Bsv^2*x^2 + Bsv*x^3)*av^4 + 2*Chi0v*x*(2*Bsv*x - 2*Bsv^2)*av^2 + 2*Chi0v*x)/(pi*(Bsv^2*av^2 + 1)^2*(Bsv^2*av^2 - 2*Bsv*av^2*x + av^2*x^2 + 1)^2);
%                     dmm(2,2) = 0;
%                     dmm(2,3) = Bsv^2/(pi*(Bsv^2*av^2 + 1)) - (Bsv - x)^2/(pi*(av^2*(Bsv - x)^2 + 1));
%                     dmm(3,3) = (2*Chi0v*av*(Bsv - x)^4)/(pi*(av^2*(Bsv - x)^2 + 1)^2) - (2*Bsv^4*Chi0v*av)/(pi*(Bsv^2*av^2 + 1)^2);

                    for i = 2:3
                        for j = (i+1):3
                            dmm(i,j) = dmm(j,i);
                        end
                    end
                    
%                     for i = 1:3
%                         g(i) = g(i) + ek*dm(i)/H(k)^2;
%                         for j = 1:3
%                             h(i,j) = h(i,j) + (dm(i)*dm(j)+ek*dmm(i,j))/H(k)^2;
%                         end
%                     end

                    for i = 1:3
                        g(i) = g(i) + ek*dm(i)/(m*log(H(k))^2);
                        for j = 1:3
                            h(i,j) = h(i,j) + (dmm(i,j)*ek/m - ek/m^2*dm(i)*dm(j)+dm(i)*dm(j)/m^2) / log(H(k))^2;
                        end
                    end
                end
                
                delta = h\g;
                delta = 0.1*delta;
                Bsv    = Bsv - delta(1);
                Chi0v  = Chi0v - delta(2);
                av     = av - delta(3);
                
                tol = sqrt((delta(1)/Bsv)^2 + (delta(2)/Chi0v)^2 + (delta(3)/av)^2);
            end
            
            assert(Chi0v*mu_o < 1);
            
            this.Bs   = Bsv;
            this.Chi0 = Chi0v;
            this.a    = av;
            
            this.CoreLossCoefficients = this.calculateCoreLossCoefficients(this.CoreLossData);
        end
        
        %% Material Properties
        function s = elementSigma(this, ~)
            s = this.Conductivity;
        end
        
        function d = elementDensity(this,~)
            d = this.Density;
        end
        
        %% Nonlinear Function Methods
        function [H, dHdB] = magnitudeH(this, B)
            [H, dHdB] = magnitudeM(this, B);
            H = B / mu_o - H;
            dHdB = 1 / mu_o - dHdB;
        end
        
        function [M, dMdB] = magnitudeM(this, B)
            Bsv   = this.Bs;
            Chi0v = this.Chi0;
            av    = this.a;
            
            M    = (Chi0v.*B)./2 + (Bsv.*Chi0v.*erf(Bsv.*av))./2 - (Chi0v.*erf(av.*(Bsv - B)).*(Bsv - B))./2 + Chi0v./(2.*pi.^(1./2).*av.*exp(Bsv.^2.*av.^2)) - Chi0v./(2.*pi.^(1./2).*av.*exp(av.^2.*(Bsv - B).^2));
            dMdB = (Chi0v.*(erf(av.*(Bsv - B)) + 1))./2;
            
            %M    = Chi0v.*(B./2 - (atan(av.*(B - Bsv)).*(B - Bsv))./pi) + (Bsv.*Chi0v.*atan(Bsv.*av))./pi;
            %dMdB = -Chi0v.*(atan(av.*(B - Bsv))./pi + (av.*(B - Bsv))./(pi.*(av.^2.*(B - Bsv).^2 + 1)) - 1./2);
        end
    end
    
    methods (Static)
        function coefficients = calculateCoreLossCoefficients(coreLossData)
            p = coreLossData(1, :).';
            f = coreLossData(2, :).';
            B = coreLossData(3, :).';

            if (numel(p) == 1 && numel(f) == 1 && numel(B) == 1)
                coefficients = [p, f, B];
            elseif ~(isempty(p) || isempty(B) || isempty(f))
                A = [ones(size(B)), log(f), log(B)];
                b = log(p);

                coefficients    = (A.'*A) \ (A.'*b);
                coefficients(1) = exp(coefficients(1));
            else
                coefficients = [0 0 0];
            end
        end
    end
end