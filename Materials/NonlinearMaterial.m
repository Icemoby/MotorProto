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
        
        BasisDegree
        
        %% Empirical Loss Data
        CoreLossData
    end
    
    properties (SetAccess = protected)
        Linear = false;
        
        %% Nonlinear Data
        MData
        
        %% Spline Data
        KnotVector
        Bbins
        Hbins
        Mbins
        SplinePFC %Coefficients of BSpline in power form
        
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
            for i = 1:(length(M)-1)
                j = i + 1;
                while (j < length(M)) && (M(j) <= M(i))
                    j = j + 1;
                end

                dMdB = (M(j) - M(i)) / (B(j) - B(i));

                if dMdB <= 0
                    for k = (i+1):j
                        M(k) = M(i)*(1+(k-i)*eps^2);
                    end
                else
                    for k = (i+1):(j-1)
                        M(k) = M(i) + dMdB * (B(k) - B(i));
                    end
                end
            end

            if M(end) < M(end-1)
                M(end) = M(end-1) * (1+eps^2);
            elseif abs(M(end)-M(end-1)) > sqrt(eps) * M(end)
                M(end+1) = M(end) * (1+eps^2);
                B(end+1) = B(end) + H(end)*mu_o;
            end

            H = B/mu_o - M;
            
            this.HData = H;
            this.BData = B;
            this.MData = M;
            
            %% Create Knot-Vector
            p = this.BasisDegree;
            k = numel(B) - 1 - p;
            U = [0*ones(1,1+p), 1:k, (k+1)*ones(1,1+p)];
            
            this.KnotVector = U;
            
            %% Calculate Knot Vector Bins
            u = unique(U)';
            w = this.evaluate_basis_funs(u,U,p);

            bins       = w*H;
            bins(end)  = H(end);
            this.Hbins = bins;

            bins       = w*M;
            bins(end)  = M(end);
            this.Mbins = bins;
            
            bins       = w*B;
            bins(end)  = B(end);
            this.Bbins = bins;
            
            %% Calculate Basis Coefficients in Power Form    
            u      = linspace(0,k+1,2*k+3)';
            w      = this.evaluate_basis_funs(u,U,p);
            w      = w*B;
            w(end) = B(end);
         	pfc    = zeros(p+1,k);
            for j = 1:(k+1)
                i = 2*j-1;
                I = i:(i+2);
                pfc(:,j) = polyfit(u(I),w(I),p);
                if p == 2 && (abs(pfc(1,j)) < sqrt(eps) * abs(pfc(2,j))) %points are colinear
                    pfc(1,j)   = 0;
                    pfc(2:3,j) = polyfit(u(I),w(I),1);
                end
            end
            pfc(p+1,1) = 0;
            
            this.SplinePFC = pfc;
            
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
            if this.BasisDegree == 1
                M = zeros(size(B));
                dMdB = zeros(size(B));
                bdata = this.BData;
                mdata = this.MData;
                n = numel(bdata);
                for i = 1:(n-1)
                    m = (mdata(i+1)-mdata(i)) / (bdata(i+1) - bdata(i));
                    I = (bdata(i) <= B) & (B < bdata(i+1));
                    M(I) = m * B(I) + (mdata(i) - m * bdata(i));
                    dMdB(I) = m;
                end
                I = (bdata(n) <= B);
                M(I) = mdata(n);
            else
                u       = this.calculate_parameter(B, this.Bbins, this.SplinePFC, this.BasisDegree);
                [w, dw] = this.evaluate_basis_funs(u, this.KnotVector, this.BasisDegree);
                M       = w * this.MData;
                dMdB    = (dw * this.MData) ./ (dw * this.BData);
            end
        end
    end
    
    methods (Static)
        function [w, dw] = evaluate_basis_funs(u, U, P)
            n = numel(U) - 1 - P;
            m = numel(u);
            
            w  = zeros(m, n);
            dw = zeros(m, n);
            for i = (1+P):n
                I = (U(i) <= u) & (u < U(i+1));
                w(I,i) = 1.0;
            end
            
            for p = 1:(P-1)
                N = w;

                i = P+1-p;
                w(:,i) = (U(i+p+1) - u) .* N(:,i+1) / (U(i+p+1) - U(i+1));

                for i = (P+2-p):(n-1)
                    w(:,i) = (u - U(i)) .* N(:,i) / (U(i+p) - U(i)) + (U(i+p+1) - u) .* N(:,i+1) / (U(i+p+1) - U(i+1));
                end

                i = n;
                w(:,i) = (u - U(i)) .* N(:,i) / (U(i+p) - U(i));
            end
            
            p = P;
            N = w;
            
            i = 1;
            w(:,i) = (U(i+p+1) - u) .* N(:,i+1) / (U(i+p+1) - U(i+1));
            for i = 2:(n-1)
                w(:,i) = (u - U(i)) .* N(:,i) / (U(i+p) - U(i)) + (U(i+p+1) - u) .* N(:,i+1) / (U(i+p+1) - U(i+1));
                dw(:,i) = p * N(:,i) / (U(i+p) - U(i)) - p * N(:,i+1) / (U(i+p+1) - U(i+1));
            end
            i = n;
            w(:,i) = (u - U(i)) .* N(:,i) / (U(i+p) - U(i));
            dw(:,i) = p * N(:,i) / (U(i+p) - U(i));


            I = (u >= U(end));
            w(I,n) = 1;
            dw(I,n) = P;
            dw(I,n-1) = -P;
        end
        
        function u = calculate_parameter(v, bins, pfc, P)
            switch P
                case 1
                    u = zeros(size(v));
                    n = numel(bins);
                    for i = 1:(n-1)
                        m = pfc(1,i);
                        b = pfc(2,i);
                        
                        I = (bins(i) <= v) & (v < bins(i+1));
                        u(I) = (v(I) - b) / m;
                    end
                    I = (bins(n) <= v);
                    u(I) = n;
                case 2
                    u = zeros(size(v));
                    n = numel(bins);
                    for i = 1:(n-1)
                        a = pfc(1,i);
                        b = pfc(2,i);
                        c = pfc(3,i);

                        I = (bins(i) <= v) & (v < bins(i+1));
                        if a == 0
                            u(I) = (v(I) - c) / b;
                        else
                            u(I) = (-b + sqrt(b^2 - 4*a*(c-v(I)))) / (2*a);
                        end
                    end
                    I = (bins(n) <= v);
                    u(I) = n;
                otherwise
                    error('#TODO - Add Newton method for evaluating parameters for p > 2 (expensive)');
            end
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