classdef (Sealed) Arc < Geometry1D
    %Arc.m Creates an object representing a 1-dimensional arc
    %   G = Arc('PropertyName',propertyvalue,...) creates an object representing
    %   a 1-dimensional arc.
    %
    %   Example: Create a parameterized arc and change its properties.
    %       P = PARAMETER_LIST;
    %       P.new('rad',1);
    %       P.new('ang',pi/3);
    %       P.new('rot',-pi/2);
    %       P.new('pos',[0.5 0.5]);
    %       G1 = Arc('Radius',   'rad', 'Angle',    'ang',...
    %                'Position', 'pos', 'Rotation', 'rot');
    %       figure;subplot(1,2,1);axis equal;
    %       G1.plot;
    %
    %       P.edit('rad',0.5);
    %       P.edit('rot',pi/2);
    %       P.edit('ang',pi);
    %       P.edit('pos',[-1 -1]);
    %       G1 = rebuild(G1);
    %       subplot(1,2,2);axis equal;
    %       G1.plot;
    %
    % Arc properties:
    %   Radius - The radius of the arc
    %   Angle  - The angle subtended by the arc
    %
    % Arc methods:
    %   Arc defines no new methods.
    %
    % Arc inheritance:
    %   Arc inherts methods and properties directly from Geometry1D. See the
    %   help for Geometry1D for more information.
    %
    %   See also Geometry1D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    properties (Dependent)
        Order        
        Radius
        Angle
    end
    
    properties (SetAccess = protected)
        Base	= 'center'
    end
    
    properties (SetAccess = private)
        pRadius
        vRadius = [0 0];
        pAngle
        vAngle  = 0;
    end
    
    methods
        function this = Arc(varargin)
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
        
        function this = set.Radius(this,radiusIn)
            this.pRadius  = radiusIn;
            this.vRadius  = this.pRadius;
            if this.fCachingEnabled;
                this          = updateCache(this);
            end
        end
        
        function radius = get.Radius(this)
            radius = this.pRadius;
        end
        
        function this = set.Angle(this,angleIn)
            this.pAngle = angleIn;
            this.vAngle = this.pAngle;
            if this.fCachingEnabled;
                this          = updateCache(this);
            end
        end
        
        function radius = get.Angle(this)
            radius = this.pAngle;
        end
        
        function orderOut = get.Order(~)
            orderOut = 2;
        end
        
        function build(this,typeIn)
        	radius_ = this.vRadius;
            
            switch this.Orientation
                case true
                    angle_ = this.vAngle;
                    rotation_ = this.Rotation;
                case false
                    angle_ = -this.vAngle;
                    rotation_ = this.Rotation-angle_;
            end
            
            position_ = this.Position;
            
            switch typeIn
                case 'x'
                    this.x       = @(s)(position_(1)+radius_*cos(rotation_+angle_*s));
                    this.fx      = false;
                case 'y'
                    this.y     	= @(s)(position_(2)+radius_*sin(rotation_+angle_*s));
                    this.fy     	= false;
                case 'dx'
                    this.dx      = @(s)(-angle_*radius_*sin(rotation_+angle_*s));
                    this.fdx     = false;
                case 'dy'
                    this.dy      = @(s)(angle_*radius_*cos(rotation_+angle_*s));
                    this.fdy     = false;
                case 'dx2'
                    this.dx2     = @(s)(-angle_^2*radius_*cos(rotation_+angle_*s));
                    this.fdx2    = false;
                case 'dy2'
                    this.dy2     = @(s)(-angle_^2*radius_*sin(rotation_+angle_*s));
                    this.fdy2    = false;
            end
        end
        
        function this = refresh(this)
            warning('Arc:refresh','Replace "refresh" with "rebuild."');
            this = rebuild(this);
        end
        
        function this = rebuild(this)
            this        = this.rebuild@Geometry1D;
            this.Radius = rebuild(this.Radius);
            this.Angle  = rebuild(this.Angle);
        end
        
        function xOut   = x(this,s)
            radius = this.vRadius;
            
            switch this.Orientation
                case true
                    angle    = this.vAngle;
                    rotation = this.vRotation;
                case false
                    angle    = -this.vAngle;
                    rotation = this.vRotation-angle;
            end
            
            position = this.vPosition;
            xOut = position(1)+radius*cos(rotation+angle*s);
        end
        
        function yOut   = y(this,s)
            radius = this.vRadius;
            
            switch this.Orientation
                case true
                    angle    = this.vAngle;
                    rotation = this.vRotation;
                case false
                    angle    = -this.vAngle;
                    rotation = this.vRotation-angle;
            end
            
            position = this.vPosition;
            yOut = position(2)+radius*sin(rotation+angle*s);
        end
        
        function dxOut	= dx(this,s)
            radius = this.vRadius;
            
            switch this.Orientation
                case true
                    angle = this.vAngle;
                    rotation = this.vRotation;
                case false
                    angle = -this.vAngle;
                    rotation = this.vRotation-angle;
            end
            
            dxOut = -angle*radius*sin(rotation + angle*s);
        end
        
        function dyOut	= dy(this,s)
            radius = this.vRadius;
            
            switch this.Orientation
                case true
                    angle = this.vAngle;
                    rotation = this.vRotation;
                case false
                    angle = -this.vAngle;
                    rotation = this.vRotation-angle;
            end
            
            dyOut = angle*radius*cos(rotation + angle*s);
        end
        
        function lengthOut = elementLength(this)
            radius    = this.vRadius;
            angle     = this.vAngle;
            lengthOut = abs(radius*angle);
        end
        
        function nOut = elementMinEdgeNumber(this)
            nOut = ceil(14 * abs(this.vAngle) / (2 * pi));
        end
        
        function dx2Out	= dx2(this,s)
            radius_ = this.vRadius;
            
            switch this.Orientation
                case true
                    angle_ = this.vAngle;
                    rotation_ = this.Rotation;
                case false
                    angle_ = -this.vAngle;
                    rotation_ = this.Rotation-angle_;
            end
            
            dx2Out = -angle_^2*radius_*cos(rotation_+angle_*s);
        end
        
        function dy2Out	= dy2(this,s)
            radius_ = this.vRadius;
            
            switch this.Orientation
                case true
                    angle_ = this.vAngle;
                    rotation_ = this.Rotation;
                case false
                    angle_ = -this.vAngle;
                    rotation_ = this.Rotation-angle_;
            end
            
            dy2Out = -angle_^2*radius_*sin(rotation_+angle_*s);
        end
        
        function areaOut = area(this)
            angle 	= this.vAngle;
            radius 	= this.vRadius;
            areaOut = radius^2 * (angle - sin(angle))/2;
            if ~this.Orientation
                areaOut = -areaOut;
            end
        end
        
        function arcsOut = split(this,sParams)
            [nRows,nCols] = size(sParams);
            if nRows > 1
                sParams = sParams.';
                nCols  = nRows;
            end
            [nRows,nCols] = size(this);
            if nRows > 1
                this = this.';
                nCols  = nRows;
            end
            
            %% Only split where necessary
            toCopy  = cellfun('isempty',sParams);
            toSplit = ~toCopy;
            
            if any(toSplit)
                thisCopy  = this(toCopy);
                thisSplit = this(toSplit);
                sParams   = sParams(toSplit);

                %% Get Original Properties
                radii    = [thisSplit.vRadius];
                angle    = [thisSplit.vAngle];
                rotation = [thisSplit.vRotation];
                position = [thisSplit.vPosition];
                position = reshape(position,2,[]).';
                base     = {thisSplit.Base};
                orient   = [thisSplit.Orientation];

                %% Calculate Parameters
                sParams  = cellfun(@(x)(horzcat(0,x)),sParams,...
                                                        'UniformOutput',false);
                dsParams = cellfun(@(x)(diff(horzcat(x,1))),sParams,...
                                                        'UniformOutput',false);

                sDims    = cellfun('length',sParams);
                nNewArcs = sum(sDims);

                sParams  = cell2mat(sParams);
                dsParams = cell2mat(dsParams);

                %%Radius
                newRadii = arrayfun(@(x,y)(repmat(x,1,y)),radii,sDims,...
                                                     	'UniformOutput',false);
                newRadii = cell2mat(newRadii);

                %%Orientation
                newOrientation = arrayfun(@(x,y)(repmat(x,1,y)),orient,sDims,...
                                                        'UniformOutput',false);
                newOrientation = cell2mat(newOrientation);

                %%Position
                newPosition      = zeros(nNewArcs,2);
                newPosition(:,1) = cell2mat(arrayfun(@(x,y)(repmat(x,y,1)),...
                                                 	position(:,1),sDims.',...
                                                    	'UniformOutput',false)...
                                              );
                newPosition(:,2) = cell2mat(arrayfun(@(x,y)(repmat(x,y,1)),...
                                                    position(:,2),sDims.',...
                                                        'UniformOutput',false)...
                                              );
                newPosition      = num2cell(newPosition,2);

                %%Base
                newBase = arrayfun(@(x,y)(repmat(x,1,y)),base,sDims,...
                                                        'UniformOutput',false);
                newBase = [newBase{:}];

                %%Rotation and Angle
                newRotation = arrayfun(@(x,y)(repmat(x,1,y)),rotation,sDims,...
                                                        'UniformOutput',false);
                newRotation = cell2mat(newRotation);

                newAngle    = arrayfun(@(x,y)(repmat(x,1,y)),angle,sDims,...
                                                        'UniformOutput',false);
                newAngle    = cell2mat(newAngle);

                newRotation = newRotation + newAngle.*(sParams);

                newAngle    = newAngle.*(dsParams);

                arcsOut     = [copy(thisCopy),Arc(nNewArcs,...
                                              	'Radius',newRadii,...
                                              	'Base',newBase,...
                                               	'Position',newPosition,...
                                              	'Angle',newAngle,...
                                               	'Rotation',newRotation,...
                                               	'Orientation',newOrientation)];
            else
                arcsOut = copy(this);
            end
        end
        
        function flagOut = coincidenceType(this,geometryIn)
            flagOut = 'none';
            while strcmp(class(geometryIn),'Rotation1D')
                geometryIn = geometryIn.BackgroundGeometry;
            end
           	isSameClass = strcmp(class(geometryIn),'Arc');
            if isSameClass
                radius1 = this.vRadius;
                radius2 = geometryIn.vRadius;
                hasSameRadius   = abs(radius1 - radius2) < sqrt(eps);
                
                if hasSameRadius
                    position1 = this.vPosition;
                    position2 = geometryIn.vPosition;
                    hasSamePosition	= norm(position1 - position2) < sqrt(eps);
                    
                    if hasSamePosition
                        rotation1 = this.Rotation;
                        rotation2 = geometryIn.Rotation;
                        angle1    = this.vAngle;
                        angle2    = geometryIn.vAngle;
                        
                        A1        = [rotation1,rotation1+angle1];
                        A2        = [rotation2,rotation2+angle2];
                        
                        hasSameAngle = all( abs(       A1  - A2) < sqrt(eps))...
                                      |all( abs(fliplr(A1) - A2) < sqrt(eps));
                        
                      	if hasSameAngle
                            flagOut = 'coincident';
                        else
                            flagOut = 'parallel';
                        end
                    end
                end
            end
        end
         
        function boolOut = isCounterClockwise(this)
            ang     = this.vAngle;
            boolOut =  ((sign(ang) == 1) & this.Orientation)...
                   	  |((sign(ang) == -1) &~this.Orientation);
        end
        
        function [W,On,Nx,Ny,S] = inOn(this,xIn,yIn)
            %% Get properties
            radius   = this.vRadius;
            angle    = this.vAngle;
            position = this.vPosition;
            rotation = this.vRotation;
            
            %% Express points in the complex plane
            zCenter = position(1) + 1i*position(2);
            z0      = this.vX0 + 1i*this.vY0;
            z1      = this.vX1 + 1i*this.vY1;
            zMid    = (z0+z1) / 2;
            zIn     = xIn + 1i*yIn;
            
            %% Adjust for orientation
            if isCounterClockwise(this);
                deltaZ       = z0-z1;
                poleIntegral = 2*pi*1i;
            else
                deltaZ       = z1-z0;
                poleIntegral = -2*pi*1i;
            end
            
            %% Calculate vector normal to the line between the arc endpoints
            zTangent = deltaZ / abs(deltaZ); 
          	zNormal  = 1i * zTangent;
            
            %% Find points that are in the closed arc
            inCircle = abs(zCenter-zIn) < radius;
            if abs(deltaZ) < sqrt(eps)
                inClosure = true(size(zIn));
            else
                dZ             = zIn - zMid + eps * zNormal;
              	zNInnerProduct = dZ * conj(zNormal);
               	inClosure      =  real(zNInnerProduct) > sqrt(eps);
            end
            
            containsPole = inClosure & inCircle;
            %% Perform Integration
            W               = log((z1-zIn)./(z0-zIn));
            W(containsPole) = W(containsPole) + poleIntegral;
            W               = W / (2 * pi * 1i);
            %W(onClosure)    = 0.5;
            
            %% Shift inputs to origin
            xShifted = xIn - position(1);
            yShifted = yIn - position(2);
            
            radiusIn = hypot(xShifted,yShifted);
            angleIn  = atan2(yShifted,xShifted);
            
            %% determine if the points lie within the angle of the arc
            if angle > 0
                if rotation+angle > 0
                    needsMod          = angleIn < rotation - sqrt(eps);
                    angleIn(needsMod) = angleIn(needsMod) + 2*pi;
                    onAngle           = angleIn < (rotation + angle + sqrt(eps));
                else
                    needsMod          = angleIn > 0;
                    angleIn(needsMod) = angleIn(needsMod) - 2*pi;
                    onAngle           = angleIn > rotation - sqrt(eps) &...
                                        angleIn < rotation + angle + sqrt(eps);
                end
            elseif angle < 0
                if rotation+angle < 0
                    needsMod          = angleIn > rotation + sqrt(eps);
                    angleIn(needsMod) = angleIn(needsMod) - 2*pi;
                    onAngle           = angleIn > (rotation + angle - sqrt(eps));
                else
                    needsMod          = angleIn < 0;
                    angleIn(needsMod) = angleIn(needsMod) + 2*pi;
                    onAngle           = angleIn < rotation + sqrt(eps) &...
                                        angleIn > rotation + angle - sqrt(eps);
                end
            else
                error('Arc:inOn','Arc angle must not be zero');
            end
            
            %% determine if the points are inside or outside the arc radius
            if isCounterClockwise(this)
                onRadius = abs(radiusIn - radius) < sqrt(eps);
                On       = onRadius & onAngle;
                Nx       = -cos(angleIn);
                Ny       = -sin(angleIn);
            else
                onRadius = abs(radiusIn - radius) < sqrt(eps);
                On       = onRadius & onAngle;
                Nx       = cos(angleIn);
                Ny       = sin(angleIn);
            end
            
          	Nx(~On) = 0;
        	Ny(~On) = 0;
            
            %% calculate S parameter
            S           = (angleIn - rotation) / angle;
            
            %% snap to end points to avoid numerical difficulties
            nearZero    = abs(S)   < sqrt(eps);
            nearOne     = abs(S-1) < sqrt(eps);
            
            S(nearZero) = 0;
            S(nearOne)  = 1;
        end
        
        function [S,On]	= cart2s(this,xIn,yIn)
            [S,On] = this.getParameter(1,xIn,yIn);
        end
        
        function [S,On] = getParameter(this,I,xIn,yIn)
            %% Get Data
            nThis = numel(this);
            rad   = [this.vRadius].';
            pos   = [this.vPosition];
            pos   = [pos(1:2:(2*nThis)); pos(2:2:(2*nThis))].';
            ang   = [this.vAngle].';
            rot   = [this.vRotation].';
            
            rad   = rad(I,:);
            pos   = pos(I,:);
            ang   = ang(I,:);
            rot   = rot(I,:);
            
         	%% Shift and rotate inputs
            xShifted = bsxfun(@minus,xIn,pos(:,1));
            yShifted = bsxfun(@minus,yIn,pos(:,2));

            radiusIn = hypot(xShifted,yShifted);
            angleIn  = atan2(yShifted,xShifted);

            minAngle = min(rot,rot+ang);
            maxAngle = max(rot,rot+ang);
            angleIn  = bsxfun(@minus,angleIn,minAngle);
            angleIn  = angleIn + sqrt(eps)*2*pi;
            angleIn  = mod(angleIn,2*pi);
            angleIn  = bsxfun(@plus,angleIn,minAngle);
            angleIn  = angleIn - sqrt(eps)*2*pi;

            On       = bsxfun(@lt,angleIn,maxAngle + sqrt(eps) * 2*pi);
            dR       = bsxfun(@minus,radiusIn,rad);
            onRadius = abs(dR) < sqrt(eps) * max(max(radiusIn));

            S        = bsxfun(@minus,angleIn,rot);
            S        = bsxfun(@rdivide,S,ang);
            On       = On & onRadius;

            %% snap to end points to avoid numerical difficulties
            nearZero    = abs(S)   < sqrt(eps);
            nearOne     = abs(S-1) < sqrt(eps);

            S(nearZero) = 0;
            S(nearOne)  = 1;

            On          = On & S >= 0 & S <= 1;
        end
    end
    
    methods (Access = protected)
        function this = updateCache(this)
            nThis =  numel(this);
            for iThis = 1:nThis
                xCache = this(iThis).x([0 0.5 1]);
                yCache = this(iThis).y([0 0.5 1]);

                this(iThis).vX0  = xCache(1);
                this(iThis).XMid = xCache(2);
                this(iThis).vX1  = xCache(3);

                this(iThis).vY0  = yCache(1);
                this(iThis).YMid = yCache(2);
                this(iThis).vY1  = yCache(3);

                if abs(this(iThis).vAngle) > pi*(1-sqrt(eps))
                    this(iThis).bbRadius = this(iThis).vRadius;
                    this(iThis).bbCenter = this(iThis).vPosition;
                else
                    this(iThis).bbRadius = sqrt( (xCache(3)-xCache(1))^2 ...
                                         +(yCache(3)-yCache(1))^2);

                    this(iThis).bbCenter = [xCache(1)+xCache(3),...
                                            yCache(1)+yCache(3)] / 2;
                end
            end
        end
    end
end
