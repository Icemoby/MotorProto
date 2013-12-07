classdef NonlinearMaterial < MaterialProperty        
    %NonlinearMaterial.m A class representing materials with nonlinear magnetic properties
    %
    % NonlinearMaterial properties:
    %   HData                    - Magnetic field intensity data from the B-H curve 
    %   BData                    - Magnetic flux density data from the B-H curve 
    %   MData                    - Magnetization calculated from HData and BData
    %   DefaultInterpolationType - Default interpolation method for the B-H curve
    %   InterpolationType        - Active interpolation method for the B-H curve
    %   HBCoefficients           - Coefficients for the active B-H curve InterpolationType
    %   CoreLossData             - LossDensity/Frequency/FluxDensity core-loss measurements
    %   CoreLossCoefficients     - Coefficients for Steinmetz equation fit to the CoreLossData
    %
    % NonlinearMaterial methods:
    %   calculateHBCoefficients            - Calculates the HBCoefficients based on the InterpolationType
    %   linearInterpolationCoefficients    -
    %   quadraticInterpolationCoefficients -
    %   analyticInterpolationCoefficients  -
    %   calculateCoreLossCoefficients      - Calculates the CoreLossCoefficients based on the CoreLossData
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
    
    %DefaultInterpolationType - 
    %
    % See also NonlinearMaterial
    DefaultInterpolationType;
    
    %InterpolationType - 
    %
    % See also NonlinearMaterial
    InterpolationType;
    
    %HBCoefficients - 
    %
    % See also NonlinearMaterial
    HBCoefficients;
    
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
        
        DefaultInterpolationType
        
        %% Empirical Loss Data
        CoreLossData
    end
    
    properties (Dependent)
        InterpolationType
    end
    
    properties (Access = protected, Hidden)
        InterpolationType_ = [];
    end

    properties (SetAccess = protected)
        Linear = false;
        
        %% Nonlinear Data
        MData
        HBCoefficients
        CoreLossCoefficients
    end
    
    methods 
        %% Constructor Methods
        function this = NonlinearMaterial(varargin)
            this = this@MaterialProperty(varargin{:});
            
            if isempty(this.InterpolationType)
                this.InterpolationType = this.DefaultInterpolationType;
            end
            this.MData                = this.calculateMData(this.HData, this.BData);
            this.CoreLossCoefficients = this.calculateCoreLossCoefficients(this.CoreLossData);
        end
        
        %% Setters
        function this = set.InterpolationType(this, value)
            this.InterpolationType_ = value;
            this.HBCoefficients     = this.calculateHBCoefficients(this.HData, this.BData, this.InterpolationType_);
        end
        
        %% Getters
        function value = get.InterpolationType(this)
            value = this.InterpolationType_;
        end
        
        %% Material Properties
        function s = elementSigma(this, ~)
            s = this.Conductivity;
        end
        
        %% Nonlinear Function Methods
        function [B, dBdH] = magnitudeB(this, H)
            switch lower(this.InterpolationType)
                case 'linear'
                    bData   = reshape(this.HData, 1, 1, []);
                    coef    = this.HBCoefficients;
                    nRegion = length(bData);

                    isGreater = bsxfun(@ge, H, bData);
                    isLess    = bsxfun(@lt, H, bData);
                    inRegion  = cat(3, isGreater(:,:,1:(end-1)) & isLess(:,:,2:end), isGreater(:,:,end));

                    slope     = zeros(size(H));
                    intercept = zeros(size(H));

                    for iRegion = 1:nRegion
                        I            = inRegion(:, :, iRegion);
                        slope(I)     = coef(iRegion, 1);
                        intercept(I) = coef(iRegion, 2);
                    end

                    B    = (H + intercept) ./ (1/mu_o - slope);
                    dBdH = 1 ./ (1/mu_o - slope);
                case 'quadratic'
                    hData     = this.HBData(1, :);
                    coef      = this.HBCoefficients;
                    nRegion   = length(hData);

                    isGreater = bsxfun(@ge, H, hData);
                    isLess    = bsxfun(@lt, H, hData);
                    inRegion  = [isGreater(:, 1:(end-1)) & isLess(:, 2:end), isGreater(:, end)];

                    B         = zeros(size(H));
                    dBdH      = zeros(size(H));
                    for iRegion = 1:nRegion
                        I       = inRegion(:, iRegion);
                        B(I)    = coef(iRegion,1)*H(I).^2 + coef(iRegion,2)*H(I) + coef(iRegion,3);
                        dBdH(I) = 2 * coef(iRegion,1)*H(I) + coef(iRegion,2);
                    end
                case 'analytic'
                    error('MotorProto:NonlinearMaterial', 'No explicit formula given');
            end
        end
        
        function [H, dHdB] = magnitudeH(this, B)
            switch lower(this.InterpolationType)
                case 'linear'
                    bData   = reshape(this.BData, 1, 1, []);
                    coef    = this.HBCoefficients;
                    nRegion = length(bData);

                    isGreater = bsxfun(@ge, B, bData);
                    isLess    = bsxfun(@lt, B, bData);
                    inRegion  = cat(3, isGreater(:,:,1:(end-1)) & isLess(:,:,2:end), isGreater(:,:,end));

                    slope     = zeros(size(B));
                    intercept = zeros(size(B));

                    for iRegion = 1:nRegion
                        I            = inRegion(:, :, iRegion);
                        slope(I)     = coef(iRegion, 1);
                        intercept(I) = coef(iRegion, 2);
                    end

                    H    = B / mu_o - slope.*B - intercept;
                    dHdB = 1 / mu_o - slope;
                case 'quadratic'
                    bData   = this.HBData(2, :);
                    coef    = this.HBCoefficients;
                    nRegion = length(bData);

                    isGreater = bsxfun(@ge,B,bData);
                    isLess    = bsxfun(@lt,B,bData);
                    inRegion  = [isGreater(:,1:(end-1)) & isLess(:,2:end), isGreater(:,end)];

                    H         = zeros(size(B));
                    dHdB      = zeros(size(B));

                    for iRegion = 1:(nRegion-1)
                        I       = inRegion(:,iRegion);
                        rad     = sqrt(coef(iRegion, 2)^2 - 4 * coef(iRegion, 1) * (coef(iRegion, 3) - B(I)));

                        H(I)    = (rad - coef(iRegion, 2)) / (2 * coef(iRegion, 1));
                        dHdB(I) = 1 ./ rad;
                    end
                    I       = inRegion(:, end);
                    H(I)    = (B(I) - coef(end, 3)) / coef(end, 2);
                    dHdB(I) = 1 ./ coef(end, 2);
                case 'analytic'
                    c    = this.HBCoefficients;
                    H    = B / mu_o - (-c(1)*c(3)*(log(exp(-2*(B-c(2))/c(3))+1)-c(4)*exp(-B/c(3))+c(4)*exp(-c(2)/c(3))*atan(exp(-(B-c(2))/c(3)))*2) + c(5));
                    dHdB = 1 / mu_o - (c(1) * (1 - tanh((B - c(2))/c(3)))*(1-c(4)*exp(-B/c(3))));
            end
        end
        
        function [M, dMdB] = magnitudeM(this, B)
            switch lower(this.InterpolationType)
                case 'linear'
                    bData   = reshape(this.BData, 1, 1, []);
                    coef    = this.HBCoefficients;
                    nRegion = length(bData);

                    isGreater = bsxfun(@ge,B,bData);
                    isLess    = bsxfun(@lt,B,bData);
                    inRegion  = cat(3, isGreater(:,:,1:(end-1)) & isLess(:,:,2:end),...
                                       isGreater(:,:,end));

                    slope     = zeros(size(B));
                    intercept = zeros(size(B));

                    for iRegion = 1:nRegion
                        I            = inRegion(:,:,iRegion);
                        slope(I)     = coef(iRegion,1);
                        intercept(I) = coef(iRegion,2);
                    end

                    M     = slope.*B + intercept;
                    dMdB    = slope;
                case 'quadratic'
                    bData   = reshape(this.HBData(2,:), 1, 1, []);
                    coef    = this.HBCoefficients;
                    nRegion = length(bData);

                    isGreater = bsxfun(@ge, B, bData);
                    isLess    = bsxfun(@lt, B, bData);
                    inRegion  = cat(3, isGreater(:, :, 1:(end-1)) & isLess(:, :, 2:end), isGreater(:, :, end));

                    H         = zeros(size(B));
                    dHdB      = zeros(size(B));

                    for iRegion = 1:(nRegion-1)
                        I       = inRegion(:, :, iRegion);
                        rad     = sqrt(coef(iRegion, 2) ^ 2 - 4 * coef(iRegion, 1) * (coef(iRegion, 3) - B(I)));

                        H(I)    = (rad - coef(iRegion, 2)) / (2 * coef(iRegion, 1));
                        dHdB(I) = 1 ./ rad;
                    end
                    I           = inRegion(:, :, end);
                    H(I)        = (B(I) - coef(end, 3)) / coef(end, 2);
                    dHdB(I)     = 1 ./ coef(end, 2);

                    M           = B / mu_o - H;
                    dMdB        = 1 / mu_o - dHdB;
                case 'analytic'
                    c    = this.HBCoefficients;
                    M    = -c(1)*c(3)*(log(exp(-2*(B-c(2))/c(3))+1)-c(4)*exp(-B/c(3))+c(4)*exp(-c(2)/c(3))*atan(exp(-(B-c(2))/c(3)))*2) + c(5);
                    I    = isinf(M);
                    M(I) = max(M(~I));
                    dMdB = c(1) * (1 - tanh((B - c(2))/c(3))).*(1-c(4)*exp(-B/c(3)));
            end
        end
    end
    
    methods (Static)
        function coefficients = calculateHBCoefficients(hData, bData, intType)
            mData = bData / mu_o - hData;
            
            switch lower(intType)
                case 'linear'
                    coefficients = NonlinearMaterial.linearInterpolationCoefficients(bData, mData);
                case 'quadratic'
                    coefficients = NonlinearMaterial.quadraticInterpolationCoefficients(hData, bData);
                case 'analytic'
                    coefficients = NonlinearMaterial.analyticInterpolationCoefficients(bData, mData);
            end
        end
        
        function coefficients = linearInterpolationCoefficients(bData, mData)
            if isrow(bData)
                bData = bData.';
            end
            
            if isrow(mData)
                mData = mData.';
            end
            
            nCoef                      = length(bData);
            coefficients               = zeros(nCoef, 2);
            coefficients(1:(end-1), 1) = (mData(2:end) - mData(1:(end-1))) ./ (bData(2:end) - bData(1:(end-1)));
            coefficients(end, 1)       = 0;
            coefficients(:, 2)         = mData - coefficients(:, 1).*bData;
        end

        function coefficients = quadraticInterpolationCoefficients(hData, bData)
            N                   = length(bData);
            coefficients        = zeros(N,3);
            coefficients(end,3) = bData(end)-hData(end)*mu_o;
            coefficients(end,2) = mu_o;
            
            for i = (N-1):-1:1
                b = [ bData(i+1);
                      bData(i);
                      2*coefficients(i+1,1)*hData(i+1) + coefficients(i+1,2)];

                A = [ hData(i+1)^2 hData(i+1) 1;
                      hData(i)^2   hData(i)   1;
                      2*hData(i+1) 1          0];

                S = sparse(1:3, 1:3, 1./sqrt(sum(abs(A),2)));
                
                A = S*A*S;
                
                coefficients(i, :) = (S*(A\(S*b))).';
            end
        end
        
        function coefOut      = analyticInterpolationCoefficients(bData, mData)
            %% Fits Magnetization data to an analytic curve whose derivative is related to tanh(x)
            %   dMdB(B) = p1 * (1 - tanh((B - p2)/p3)) * (1 - p4*exp(-B/p3))

            M = mData(2:end);
            B = bData(2:end);
            M = reshape(M,1,1,[]);
            B = reshape(B,1,1,[]);
            
            %% Minimize Over A Discrete Set Of Parameters To Get A Good Initial Condition
            ND = 11;
            p1 = linspace(min(M./B), max(M./B), ND);
            p2 = linspace(1.5,       2,         ND);
            p3 = linspace(0,         0.5,       ND);
            p4 = linspace(0,         1e-3,      ND); 
            
            f  = zeros(numel(p1), numel(p2), numel(p3), numel(p4));
            for i = 1:numel(p1)
                for j = 1:numel(p2)
                    for k = 1:numel(p3)
                        for l = 1:numel(p4)
                            f(i,j,k,l) = sum(1./M.^2.*(M + (p1(i).*p3(k) .*( p4(l) + log(exp(-(B.*2-p2(j).*2)./p3(k))+1) - log(exp((p2(j).*2)./p3(k))+1) - p4(l).*(  exp(-B./p3(k)) + exp(-p2(j)./p3(k)).*atan(exp(p2(j)./p3(k))).*2 - exp(-p2(j)./p3(k)).*atan(exp(-(B-p2(j))./p3(k))).*2)))./(tanh(p2(j)./p3(k))+1)).^2,...
                                             3);
                        end
                    end
                end
            end
            
            [f,l] = min(f,[],4);
            [f,k] = min(f,[],3);
            [f,j] = min(f,[],2);
            [f,i] = min(f,[],1);
            
            j  = j(i);
            k  = k(i,j);
            l  = l(i,j,k);
            
            p1 = p1(i);
            p2 = p2(j);
            p3 = p3(k);
            p4 = p4(l);
            
            %% Use Newton's Method to find the local minimum
            maxIter = 100;
            fOut  	= zeros(maxIter, 1);
            M       = reshape(M, 1, 1, []);
            B       = reshape(B, 1, 1, []);
            dF      = inf;
            fOut(1) = f;
            gOut    = ones(3,1);
            iIter   = 1;
            while iIter < maxIter && dF > sqrt(eps) && norm(gOut) > sqrt(eps)
                gOut = sum([(1.0./M.^2.*p3.*(M+(p1.*p3.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)).*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0).*2.0)./(tanh(p2./p3)+1.0);1.0./M.^2.*(M+(p1.*p3.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)).*(p1.*1.0./(tanh(p2./p3)+1.0).^2.*(tanh(p2./p3).^2-1.0).*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0)-(p1.*p3.*((p4.*2.0)./(p3.*(exp((p2.*2.0)./p3)+1.0))-(exp(-(B.*2.0-p2.*2.0)./p3).*2.0)./(p3.*(exp(-(B.*2.0-p2.*2.0)./p3)+1.0))+(exp((p2.*2.0)./p3).*2.0)./(p3.*(exp((p2.*2.0)./p3)+1.0))-(p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0)./p3+(p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0)./p3-(p4.*exp(-p2./p3).*exp(-(B-p2)./p3).*2.0)./(p3.*(exp((B.*-2.0+p2.*2.0)./p3)+1.0))))./(tanh(p2./p3)+1.0)).*2.0;1.0./M.^2.*(M+(p1.*p3.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)).*((p1.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)+(p1.*p3.*((p2.*1.0./p3.^2.*exp((p2.*2.0)./p3).*2.0)./(exp((p2.*2.0)./p3)+1.0)+(1.0./p3.^2.*exp(-(B.*2.0-p2.*2.0)./p3).*(B.*2.0-p2.*2.0))./(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-B.*1.0./p3.^2.*p4.*exp(-B./p3)+(p2.*1.0./p3.^2.*p4.*2.0)./(exp((p2.*2.0)./p3)+1.0)-p2.*1.0./p3.^2.*p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p2.*1.0./p3.^2.*p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0+(1.0./p3.^2.*p4.*exp(-p2./p3).*exp(-(B-p2)./p3).*(B-p2).*2.0)./(exp((B.*-2.0+p2.*2.0)./p3)+1.0)))./(tanh(p2./p3)+1.0)-(p1.*p2.*1.0./(tanh(p2./p3)+1.0).^2.*(tanh(p2./p3).^2-1.0).*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./p3).*2.0;(1.0./M.^2.*p1.*p3.*(M+(p1.*p3.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)).*(exp(-B./p3)+exp(-p2./p3).*atan(exp(p2./p3)).*2.0-exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0-1.0).*-2.0)./(tanh(p2./p3)+1.0)],3);
                
                hOut = sum([[(2.*p3.^2.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)).^2)./(M.^2.*(tanh(p2./p3) + 1).^2), (2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1).^2) - (2.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(M.^2.*(tanh(p2./p3) + 1)) + (2.*p3.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1)), (2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1)) + (2.*p3.*((p1.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1)) + (2.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*p3.*(tanh(p2./p3) + 1).^2), - (2.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*p3.^2.*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1).^2)];
                            [(2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(((tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)))./M.^2 + (2.*p3.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1)), (2.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).^2)./M.^2 - (2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((2.*p1.*(tanh(p2./p3).^2 - 1).*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) + 4./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - 4./(p3.^2.*exp((2.*(2.*B - 2.*p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2) - (4.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (4.*exp((4.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1).^2) - (2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (4.*p4.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1).^2) - (2.*p4)./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1)) - (4.*p4)./(p3.^2.*exp(p2./p3).*exp((3.*(B - p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2)))./(tanh(p2./p3) + 1) - (2.*p1.*(tanh(p2./p3).^2 - 1).^2.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^3) + (2.*p1.*tanh(p2./p3).*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)))./M.^2, (2.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*((p1.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)))./M.^2 - (2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p1.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*(tanh(p2./p3).^2 - 1).*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1).^2 + (p1.*p3.*(2./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) - (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) - (4.*p2.*exp((2.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1)) + (4.*p2.*exp((4.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1).^2) - (2.*(2.*B - 2.*p2))./(p3.^3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*(2.*B - 2.*p2))./(p3.^3.*exp((2.*(2.*B - 2.*p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2) + (2.*p2.*p4)./(p3.^3.*(exp((2.*p2)./p3) + 1)) + (2.*p4)./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^3.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^3.*exp(p2./p3)) + (4.*p2.*p4.*exp((2.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1).^2) - (2.*p2.*p4)./(p3.^3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1)) + (4.*p4.*(B - p2))./(p3.^3.*exp(p2./p3).*exp((3.*(B - p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2)))./(tanh(p2./p3) + 1) + (2.*p1.*p2.*(tanh(p2./p3).^2 - 1).^2.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.^2.*(tanh(p2./p3) + 1).^3) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(p3.*(tanh(p2./p3) + 1).^2) - (2.*p1.*p2.*tanh(p2./p3).*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.^2.*(tanh(p2./p3) + 1).^2)))./M.^2, - (2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p1.*(tanh(p2./p3).^2 - 1).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(tanh(p2./p3) + 1).^2 + (p1.*p3.*(2./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - 2./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)))./M.^2 - (2.*p1.*p3.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1))];
                            [(2.*(M + (p1.*p3.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3))./(tanh(p2./p3) + 1) + (p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)))./M.^2 + (2.*p3.*((p1.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1)), (2.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*((p1.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)))./M.^2 - (2.*(M + (p1.*p3.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p1.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*(tanh(p2./p3).^2 - 1).*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1).^2 + (p1.*p3.*(2./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) - (2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) - (4.*p2.*exp((2.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1)) + (4.*p2.*exp((4.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1).^2) - (2.*(2.*B - 2.*p2))./(p3.^3.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) + (2.*(2.*B - 2.*p2))./(p3.^3.*exp((2.*(2.*B - 2.*p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2) + (2.*p2.*p4)./(p3.^3.*(exp((2.*p2)./p3) + 1)) + (2.*p4)./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^3.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^3.*exp(p2./p3)) + (4.*p2.*p4.*exp((2.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1).^2) - (2.*p2.*p4)./(p3.^3.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1)) + (4.*p4.*(B - p2))./(p3.^3.*exp(p2./p3).*exp((3.*(B - p2))./p3).*(exp(-(2.*(B - p2))./p3) + 1).^2)))./(tanh(p2./p3) + 1) + (2.*p1.*p2.*(tanh(p2./p3).^2 - 1).^2.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.^2.*(tanh(p2./p3) + 1).^3) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(p3.*(tanh(p2./p3) + 1).^2) - (2.*p1.*p2.*tanh(p2./p3).*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.^2.*(tanh(p2./p3) + 1).^2)))./M.^2, (2.*((p1.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)).^2)./M.^2 - (2.*(M + (p1.*p3.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p1.*p3.*((B.^2.*p4)./(p3.^4.*exp(B./p3)) - (2.*p2.^2.*p4)./(p3.^4.*(exp((2.*p2)./p3) + 1)) + (4.*p2.*exp((2.*p2)./p3))./(p3.^3.*(exp((2.*p2)./p3) + 1)) + (2.*(2.*B - 2.*p2))./(p3.^3.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (2.*B.*p4)./(p3.^3.*exp(B./p3)) + (4.*p2.^2.*exp((2.*p2)./p3))./(p3.^4.*(exp((2.*p2)./p3) + 1)) - (4.*p2.^2.*exp((4.*p2)./p3))./(p3.^4.*(exp((2.*p2)./p3) + 1).^2) + (4.*p2.*p4)./(p3.^3.*(exp((2.*p2)./p3) + 1)) - (2.*B - 2.*p2).^2./(p3.^4.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) + (2.*B - 2.*p2).^2./(p3.^4.*exp((2.*(2.*B - 2.*p2))./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1).^2) - (4.*p2.^2.*p4.*exp((2.*p2)./p3))./(p3.^4.*(exp((2.*p2)./p3) + 1).^2) - (4.*p2.*p4.*atan(exp(p2./p3)))./(p3.^3.*exp(p2./p3)) + (4.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^3.*exp(p2./p3)) + (2.*p2.^2.*p4.*atan(exp(p2./p3)))./(p3.^4.*exp(p2./p3)) - (2.*p2.^2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^4.*exp(p2./p3)) - (2.*p4.*(B - p2).^2)./(p3.^4.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1)) + (4.*p4.*(B - p2).^2)./(p3.^4.*exp(p2./p3).*exp((3.*(B - p2))./p3).*(exp(-(2.*(B - p2))./p3) + 1).^2) + (4.*p4.*(B - p2))./(p3.^3.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1)) - (4.*p2.*p4.*(B - p2))./(p3.^4.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (2.*p1.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) + (2.*p1.*p2.*(tanh(p2./p3).^2 - 1).*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(p3.*(tanh(p2./p3) + 1).^2) - (2.*p1.*p2.^2.*(tanh(p2./p3).^2 - 1).^2.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.^3.*(tanh(p2./p3) + 1).^3) + (2.*p1.*p2.^2.*tanh(p2./p3).*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.^3.*(tanh(p2./p3) + 1).^2)))./M.^2, (2.*(M + (p1.*p3.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((p1.*p3.*((2.*p2)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - B./(p3.^2.*exp(B./p3)) - (2.*p2.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*(exp(-B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(exp(-(B - p2)./p3)))./exp(p2./p3) - 1))./(tanh(p2./p3) + 1) + (p1.*p2.*(tanh(p2./p3).^2 - 1).*(exp(-B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(exp(-(B - p2)./p3)))./exp(p2./p3) - 1))./(p3.*(tanh(p2./p3) + 1).^2)))./M.^2 - (2.*p1.*p3.*((p1.*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(exp(-(2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(exp(-(B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(exp(-(2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(exp(-(2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(exp(-(B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)).*(exp(-B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(exp(-(B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1))];
                            [- (2.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*p3.^2.*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(M.^2.*(tanh(p2./p3) + 1).^2), - (2.*p1.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(2./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - 2./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*p3.*((p1.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1).^2 - (p1.*p3.*((2.*p4)./(p3.*(exp((2.*p2)./p3) + 1)) - 2./(p3.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) + (2.*exp((2.*p2)./p3))./(p3.*(exp((2.*p2)./p3) + 1)) - (2.*p4.*atan(exp(p2./p3)))./(p3.*exp(p2./p3)) + (2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.*exp(p2./p3)) - (2.*p4)./(p3.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(tanh(p2./p3).^2 - 1).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1).^2), (2.*p1.*p3.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*((2.*p2)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - B./(p3.^2.*exp(B./p3)) - (2.*p2.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1)) - (2.*p1.*p3.*((p1.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1) + (p1.*p3.*((2.*p2.*exp((2.*p2)./p3))./(p3.^2.*(exp((2.*p2)./p3) + 1)) + (2.*B - 2.*p2)./(p3.^2.*exp((2.*B - 2.*p2)./p3).*(1./exp((2.*B - 2.*p2)./p3) + 1)) - (B.*p4)./(p3.^2.*exp(B./p3)) + (2.*p2.*p4)./(p3.^2.*(exp((2.*p2)./p3) + 1)) - (2.*p2.*p4.*atan(exp(p2./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p2.*p4.*atan(1./exp((B - p2)./p3)))./(p3.^2.*exp(p2./p3)) + (2.*p4.*(B - p2))./(p3.^2.*exp(p2./p3).*exp((B - p2)./p3).*(1./exp((2.*(B - p2))./p3) + 1))))./(tanh(p2./p3) + 1) - (p1.*p2.*(tanh(p2./p3).^2 - 1).*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(p3.*(tanh(p2./p3) + 1).^2)).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*(tanh(p2./p3) + 1)) + (2.*p1.*p2.*(M + (p1.*p3.*(p4 + log(1./exp((2.*B - 2.*p2)./p3) + 1) - log(exp((2.*p2)./p3) + 1) - p4./exp(B./p3) - (2.*p4.*atan(exp(p2./p3)))./exp(p2./p3) + (2.*p4.*atan(1./exp((B - p2)./p3)))./exp(p2./p3)))./(tanh(p2./p3) + 1)).*(tanh(p2./p3).^2 - 1).*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1))./(M.^2.*p3.*(tanh(p2./p3) + 1).^2), (2.*p1.^2.*p3.^2.*(1./exp(B./p3) + (2.*atan(exp(p2./p3)))./exp(p2./p3) - (2.*atan(1./exp((B - p2)./p3)))./exp(p2./p3) - 1).^2)./(M.^2.*(tanh(p2./p3) + 1).^2)];]...
                            ,3);

                dP = hOut\gOut;
                
                if any(isnan(dP) | isinf(dP))
                    dP(:) = 0;
                    warning('MotorProto:NonlinearMaterial', ['Interpolation stagnated with RMS percent error %d.'...
                                                             'Values much less than 1 indicate a good fit.'         ], sqrt(fOut(iIter) / numel(M)));
                end
                
                p1 = p1 - dP(1);
                p2 = p2 - dP(2);
                p3 = p3 - dP(3);
                p4 = p4 - dP(4);
                
                if p4 < 0
                    p4 = 0;
                end
                
                iIter       = iIter+1;
                fOut(iIter) = sum(1.0./M.^2.*(M+(p1.*p3.*(p4+log(exp(-(B.*2.0-p2.*2.0)./p3)+1.0)-log(exp((p2.*2.0)./p3)+1.0)-p4.*exp(-B./p3)-p4.*exp(-p2./p3).*atan(exp(p2./p3)).*2.0+p4.*exp(-p2./p3).*atan(exp(-(B-p2)./p3)).*2.0))./(tanh(p2./p3)+1.0)).^2);
                dF          = abs(fOut(iIter)-fOut(iIter-1))/(fOut(iIter)+sqrt(eps));
            end

            p1      = p1 / (1-tanh(-p2/p3));
            p5      = -p1*p3*(p4-log(exp(p2*2/p3)+1)-p4*exp(-p2/p3)*atan(exp(p2/p3))*2);
            coefOut = [p1, p2, p3, p4, p5];
        end
        
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