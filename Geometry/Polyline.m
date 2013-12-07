classdef (Sealed) Polyline < Geometry1D
    %Polyline.m Creates an object representing a 1-dimensional line
    %   G = Polyline('PropertyName',propertyvalue,...) creates an object
    %   representing a 1-dimensional line
    %
    %   Example: Create a parameterized line and change its properties.
    %       P = PARAMETER_LIST;
    %       P.new('pts',[0 0;1 1]);
    %       P.new('pos',[0 0]);
    %       P.new('rot',0);
    %
    %       G1 = Polyline('Points','pts','Position','pos','Rotation','rot');
    %       figure;subplot(1,2,1);axis equal;
    %       G1.plot;
    %
    %       P.edit('pts',[0 1;1 0]);
    %       P.edit('pos',[1 1]);
    %       P.edit('rot',pi/8);
    %       G1 = rebuild(G1);
    %       subplot(1,2,2);axis equal;
    %       G1.plot;
    %
    % Polyline properties:
    %   Points - A 2 by 2 matrix containing endpoint coordinates
    %
    %   See also Geometry1D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    properties (Dependent)
        Order
        Points
    end
    
    properties (SetAccess = private)
        pPoints
        vPoints
    end
    
    methods
        function this = Polyline(varargin)
            this = this@Geometry1D;
            if nargin~=0
                if mod(nargin,2) == 0                    
                    for iArg = 1:2:nargin
                        this.(varargin{iArg})=varargin{iArg+1};
                    end
                else
                    thisSize         = varargin{1};
                    this(1,thisSize) = this;
                    this(:)          = this(1,end);
                    this             = copy(this);
                    for iArg = 2:2:nargin
                        if ~iscell(varargin{iArg+1})
                            varargin{iArg+1} = num2cell(varargin{iArg+1});
                        end
                        [this.(varargin{iArg})] = deal(varargin{iArg+1}{:});
                    end
                end
                this.enableCaching;
                this.updateCache;
            end
        end

        function this = set.Points(this,pointsIn)
            this.pPoints  = this.setProperty(pointsIn);
            this.vPoints  = this.pPoints.Value;
            
            if this.fCachingEnabled;
                this = updateCache(this);
            end
            
            [n,m]         = size(this.vPoints);
            assert( n==2 & m==2,'Points:Line1D',...
                    'The points property must be a 2 by 2 matrix');
        end
        
        function points = get.Points(this)
            points = this.pPoints;
        end
        
        function orderOut = get.Order(~)
            orderOut = 1;
        end
        
        function build(this,typeIn)
            %% build the anonymous function defined by param
            points_     = this.vPoints;
            if ~this.Orientation
                points_ = points_(end:-1:1,:);
            end
            
            rotation_   = this.Rotation.Value;
            rotMat      = [  cos(rotation_) sin(rotation_);
                            -sin(rotation_) cos(rotation_)];
            
            points_ = points_*rotMat;
            position_ = this.Position.Value;
            
            switch typeIn
                case 'x'
                    this.x       = @(s)(position_(1)+(1-s)*points_(1,1)+s*points_(2,1));
                    this.fx      = false;
                case 'y'
                    this.y       = @(s)(position_(2)+(1-s)*points_(1,2)+s*points_(2,2));
                    this.fy      = false;
                case 'dx'
                    this.dx      = @(s)((-points_(1,1)+points_(2,1)));
                    this.fdx     = false;
                case 'dy'
                    this.dy      = @(s)(-points_(1,2)+points_(2,2));
                    this.fdy     = false;
                case 'dx2'
                    this.dx2     = @(s)(0);
                    this.fdx2	= false;
                case 'dy2'
                    this.dy2     = @(s)(0);
                    this.fdy2    = false;
            end
        end
        
        function this = refresh(this)
            warning('Polyline:refresh','Replace "refresh" with "rebuild."');
            this = rebuild(this);
        end
        
        function this = rebuild(this)
            this          = this.rebuild@Geometry1D;
            this.Points   = rebuild(this.Points);
        end
 
    	function xOut = x(this,s)
        	points = this.vPoints;
            if ~this.Orientation
                points = flipud(points);
            end
            
            rotation       = this.vRotation;
            rotationMatrix = [  cos(rotation) sin(rotation);
                               -sin(rotation) cos(rotation)];
            
            points   = points * rotationMatrix;
            position = this.vPosition;
            xOut     = position(1) + (1-s)*points(1,1) + s*points(2,1);
        end
        
     	function yOut = y(this,s)
        	points = this.vPoints;
            if ~this.Orientation
                points = points(end:-1:1,:);
            end
            
            rotation       = this.vRotation;
            rotationMatrix = [  cos(rotation) sin(rotation);
                               -sin(rotation) cos(rotation)];
            
            points   = points*rotationMatrix;
            position = this.vPosition;
            yOut	 = position(2) + (1-s)*points(1,2) + s*points(2,2);
        end
        
        function dxOut 	= dx(this,s)
            dxOut  = zeros(size(s));
        	points = this.vPoints;
            if ~this.Orientation
                points = points(end:-1:1,:);
            end
            
            rotation       = this.vRotation;
            rotationMatrix = [  cos(rotation) sin(rotation);
                               -sin(rotation) cos(rotation)];
            
            points     = points*rotationMatrix;
            dxOut(:,:) = -points(1,1) + points(2,1);
        end
        
        function dyOut 	= dy(this,s)
            dyOut  = zeros(size(s));
        	points = this.vPoints;
            if ~this.Orientation
                points = points(end:-1:1,:);
            end
            
            rotation      = this.vRotation;
            rotationMatrix = [  cos(rotation) sin(rotation);
                               -sin(rotation) cos(rotation)];
            
            points     = points*rotationMatrix;
            dyOut(:,:) = -points(1,2) + points(2,2);
        end
        
        function lengthOut = elementLength(this)
            points    = this.vPoints;
            lengthOut = sqrt(   (points(1,1)-points(2,1)).^2 ...
                               +(points(1,2)-points(2,2)).^2);
        end
        
        function nOut = elementMinEdgeNumber(~)
            nOut = 1;
        end
        
        function dx2Out	= dx2(~,~)
            dx2Out = 0;
        end
        
        function dy2Out	= dy2(~,~)
            dy2Out = 0;
        end
        
        function areaOut = area(~)
            areaOut = 0;
        end
        
        function [S,On] = cart2s(this,xIn,yIn)
%             [~,I,~,S] = this.inOn(xIn,yIn);%
            x0 = this.vX0;
            x1 = this.vX1;
            y0 = this.vY0;
            y1 = this.vY1;
            dx = x1 - x0;
            dy = y1 - y0;
            if abs(dx)  < sqrt(eps)
                S  = (yIn - y0) / dy;
                On = abs(xIn - x0) < sqrt(eps);
            elseif abs(dy) < sqrt(eps)
                S  = (xIn - x0) / dx;
                On = abs(yIn - y0) < sqrt(eps);
            else
                Sx = (xIn - x0) / dx;
                Sy = (yIn - y0) / dy;
                S  = (Sx + Sy) / 2;
                On = abs(Sx - Sy) < sqrt(eps);
            end

            nearZero    = abs(S)   < sqrt(eps);
            nearOne     = abs(S-1) < sqrt(eps);
            S(nearZero) = 0;
            S(nearOne)  = 1;
            On          = On & S >= 0 & S <= 1;
        end
        
        function [S,On] = getParameter(this,I,xIn,yIn)    
            %% Get Data            
            X  = horzcat([this.vX0].', [this.vX1].');
            Y  = horzcat([this.vY0].', [this.vY1].');    
            
            X  = X(I,:);
            Y  = Y(I,:);
            
            dX = X(:,2) - X(:,1);
            dY = Y(:,2) - Y(:,1);
            
            %% Calculate Parameter
            [nr,nc] = size(xIn);
            isVert  = abs(dX) < sqrt(eps);
            isHorz  = abs(dY) < sqrt(eps);
            isNorm  = ~(isHorz | isVert);

            S  = zeros(nr,nc);
            On = false(nr,nc);

            Sx           = bsxfun(@minus,xIn,X(:,1));
            On(isVert,:) = abs(Sx(isVert,:)) < sqrt(eps);
            Sx           = bsxfun(@rdivide,Sx,dX);
            S(isHorz,:)  = Sx(isHorz,:);

            Sy           = bsxfun(@minus,yIn,Y(:,1));
            On(isHorz,:) = abs(Sy(isHorz,:)) < sqrt(eps);
            Sy           = bsxfun(@rdivide,Sy,dY);
            S(isVert,:)  = Sy(isVert,:);

            S(isNorm,:)  =    (Sx(isNorm,:) + Sy(isNorm,:)) / 2;
            On(isNorm,:) = abs(Sx(isNorm,:) - Sy(isNorm,:)) < sqrt(eps);
            
            %% Snap endpoints
            nearZero    = abs(S)   < sqrt(eps);
            nearOne     = abs(S-1) < sqrt(eps);
            S(nearZero) = 0;
            S(nearOne)  = 1;

            On          = On & S >= 0 & S <= 1;
        end
        
        function flagOut = coincidenceType(this,geometryIn)
            flagOut     = 'none';
            
            if isa(geometryIn,'Polyline')
                %% if slopes are the same
                points1 = [this.vX0,this.vY0;
                           this.vX1,this.vY1];
                deltax1 = points1(2,1) - points1(1,1);
                deltay1 = points1(2,2) - points1(1,2);
                
                points2 = [geometryIn.vX0,geometryIn.vY0;
                           geometryIn.vX1,geometryIn.vY1];
                deltax2 = points2(2,1) - points2(1,1);
                deltay2 = points2(2,2) - points2(1,2);
                
                slopesAreEqual = abs(deltax1*deltay2 - deltax2*deltay1)...
                                        < sqrt(eps);
                
                if slopesAreEqual
                    flagOut = 'parallel';
                    
                    %% if endpoints are equal
                    dpoints1 = max( max( abs(points1 -        points2 )) );
                    dpoints2 = max( max( abs(points1 - flipud(points2))) );
                    
                    endPointsAreSame = dpoints1 < sqrt(eps) | dpoints2 < sqrt(eps);
                    
                    if endPointsAreSame
                        flagOut = 'coincident';
                    end
                end
            end
        end
        
        function [W,On,Nx,Ny,S] = inOn(this,xIn,yIn)
            %% Get Data
            x0     = this.vX0;
            x1     = this.vX1;
            y0     = this.vY0;
            y1     = this.vY1;
            deltax = x1-x0;
            deltay = y1-y0;
            
            %% Cacluate Winding Number
            z0     = x0  + 1i*y0;
            z1     = x1  + 1i*y1;
            zk     = xIn + 1i*yIn;
            dz0    = z0 - zk;
            dz1    = z1 - zk;
            W      = log(dz1./dz0);
            W      = W / (2 * pi * 1i);
            
            %% Calculate parameter value
            if abs(deltax)  < sqrt(eps)
                S  = (yIn - y0) / deltay;
                On = abs(xIn - (x0 + x1) / 2) < sqrt(eps);
            elseif abs(deltay) < sqrt(eps)
                S  = (xIn - x0) / deltax;
                On = abs(yIn - (y0 + y1) / 2) < sqrt(eps);
            else
                S1 = (xIn - x0) / deltax;
                S2 = (yIn - y0) / deltay;
                S  = (S1 + S2) / 2;
                On = abs(S1 - S2) < sqrt(eps);
            end
            
            %% Snap to endpoints to avoid numerical difficulties
            nearZero    = abs(S)   < sqrt(eps);
            nearOne     = abs(S-1) < sqrt(eps);
            S(nearZero) = 0;
            S(nearOne)  = 1;
            
            %% Ensure points are within range
            On = On & S >= 0 & S <= 1;
            
            %% Calculate right-handed normal vector
            normalVector = [-deltay; deltax];
            normalVector = normalVector / norm(normalVector);
            
            %% Assign the normal vector to points on the line
            n       = length(xIn);
            Nx      = normalVector(1)*ones(n,1);
            Ny      = normalVector(2)*ones(n,1);
           	Nx(~On) = 0;
           	Ny(~On) = 0;
        end
        
        function linesOut = split(this,sParams)
            [nRows,~] = size(sParams);
            if nRows > 1
                sParams = sParams.';
            end
            [nRows,~] = size(this);
            if nRows > 1
                this = this.';
            end
            
            %% Only split where necessary
            toCopy  = cellfun('isempty',sParams);
            toSplit = ~toCopy;
            
            if any(toSplit)
                thisCopy  = this(toCopy);
                thisSplit = this(toSplit);
                sParams   = sParams(toSplit);
                
                %% Get Original Properties
                points   = [thisSplit.vPoints];
                rotation = [thisSplit.vRotation];
                position = [thisSplit.vPosition];
                position = reshape(position,2,[]).';
                orient   = [thisSplit.Orientation];

                X0       = points(1,1:2:end);
                X1       = points(2,1:2:end);
                Y0       = points(1,2:2:end);
                Y1       = points(2,2:2:end);

                XTemp       = X0(~orient);
                X0(~orient) = X1(~orient);
                X1(~orient) = XTemp;

                YTemp       = Y0(~orient);
                Y0(~orient) = Y1(~orient);
                Y1(~orient) = YTemp;    

                %% Calculate Parameters
                sParams   = cellfun(@(x)(horzcat(0,x)),sParams,...
                                                        'UniformOutput',false);
                sDims     = cellfun('length',sParams);
                nNewLines = sum(sDims);

                dsParams = cellfun(@(x)(diff(horzcat(x,1))),sParams,...
                                                        'UniformOutput',false);

                sParams   = cell2mat(sParams);
                dsParams  = cell2mat(dsParams);   

                %%Points
                X0 = arrayfun(@(x,y)(repmat(x,1,y)),X0,sDims,...
                                                        'UniformOutput',false);
                X0 = cell2mat(X0);

                X1 = arrayfun(@(x,y)(repmat(x,1,y)),X1,sDims,...
                                                    	'UniformOutput',false);
                X1 = cell2mat(X1);

                Y0 = arrayfun(@(x,y)(repmat(x,1,y)),Y0,sDims,...
                                                      	'UniformOutput',false);
                Y0 = cell2mat(Y0);

                Y1 = arrayfun(@(x,y)(repmat(x,1,y)),Y1,sDims,...
                                                    	'UniformOutput',false);
                Y1 = cell2mat(Y1);

                newX0 = X0 + (X1-X0).*sParams;
                newX1 = X0 + (X1-X0).*(sParams+dsParams);            

                newY0 = Y0 + (Y1-Y0).*sParams;
                newY1 = Y0 + (Y1-Y0).*(sParams+dsParams);

                newPoints = arrayfun(@(a,b,c,d)([a b;c d]),...
                                        newX0,newY0,newX1,newY1,...
                                                        'UniformOutput',false);

                %%Position
                newPosition      = zeros(nNewLines,2);
                newPosition(:,1) = cell2mat(arrayfun(@(x,y)(repmat(x,y,1)),...
                                            	position(:,1),sDims.',...
                                                        'UniformOutput',false)...
                                              );
                newPosition(:,2) = cell2mat(arrayfun(@(x,y)(repmat(x,y,1)),...
                                              	position(:,2),sDims.',...
                                                        'UniformOutput',false)...
                                              );
                newPosition      = num2cell(newPosition,2);

                %%Rotation and Angle
                newRotation = arrayfun(@(x,y)(repmat(x,1,y)),...
                                              rotation,sDims,...
                                                        'UniformOutput',false);
                newRotation = cell2mat(newRotation);

                linesOut    = [copy(thisCopy),Polyline(nNewLines,...
                                                   	'Points',newPoints,...
                                                   	'Position',newPosition,...
                                                    'Rotation',newRotation)];
            else
                linesOut = copy(this);
            end
        end
    end
    
    methods (Access = protected)
        function this = updateCache(this)
            nThis = numel(this);
            for iThis = 1:nThis
                xCache           = this(iThis).x([0 0.5 1]);
                yCache           = this(iThis).y([0 0.5 1]);

                this(iThis).vX0  = xCache(1);
                this(iThis).XMid = xCache(2);
                this(iThis).vX1  = xCache(3);

                this(iThis).vY0  = yCache(1);
                this(iThis).YMid = yCache(2);
                this(iThis).vY1  = yCache(3);

                this(iThis).bbRadius = sqrt( (xCache(3)-xCache(1))^2 ...
                                            +(yCache(3)-yCache(1))^2);

                this(iThis).bbCenter = [xCache(2),yCache(2)];
            end
        end
    end
end