classdef Intersection2D < Geometry2D
    %Intersection2D.m Calculates the intersection of two planes
    %   C = Intersection2D(G1,G2) Creates a new plane C which is the
    %   overlapping portions of G1 and G2.
    %
    %   When G1 and G2 are planes, the following are equivalent statements
    %
    %       C = Intersection2D(G1,G2)
    %       C = G1 * G2
    %       C = intersection(G1,G2)
    %       C = G1.intersection(G2)
    %
    %   The internal edges of C are automatically discarded. If G1 and G2
    %   intersect on a boundary curve or a set of boundary curves, C is an
    %   open plane containing only those boundary curves. In this case,
    %   area(C) returns the area of the convex hull defined by the
    %   boundary curves.
    %
    %   Example: Calculate the intersection of a rectangle and annullus.
    %       G1 = Geometry2D.draw('Rectangle2D',...
    %                               'Length',0.5,...
    %                               'Width',0.5,...
    %                               'Base','Corner',...
    %                               'Position',[0.5 0.5]);
    %
    %       G2 = Geometry2D.draw('Sector2D',...
    %                               'Radius',[0.5 1],...
    %                               'Angle',pi/2);
    %
    %       G3 = G1 * G2;
    %       G4 = intersection(G1,G2);
    %
    %       figure;subplot(2,2,1);axis equal;
    %       G1.plot;
    %       subplot(2,2,2);axis equal;
    %       G2.plot;
    %       subplot(2,2,3);axis equal;
    %       G3.plot;
    %       subplot(2,2,4);axis equal;
    %       G4.plot;
    %
    % Intersection2D methods:
    %   Intersection2D defines no new methods
    %
    % Intersection2D properties:
    %   Intersection2D defines no new properties
    %
    % Intersection2D inheritance:
    %   Intersection2D inherts methods and properties from Composite2D. See the 
    %   help for Composite2D for more information.
    %
    % See also Composite2D, Geometry2D, Geometry, Parameterizable, MotorProto
    
    %   Copyright 2011 Jason Pries
    %   $Revision 0.0.0.2 $
    
    properties
        InputGeometry
    end
    
    methods
        %% Constructor Methods
        function this = Intersection2D(geometry,varargin)         	
            this = this@Geometry2D;
            if nargin~=0
                gClasses        = geometry.getElementClasses;
                isIntersection  = strcmp('Intersection2D',gClasses);
                
                if any(isIntersection)
                    this.InputGeometry = [geometry(~isIntersection),...
                                          geometry(isIntersection).InputGeometry];
                else
                    this.InputGeometry = geometry;
                end

                this.PlotStyle     = geometry(1).PlotStyle;
                for i = 1:2:(nargin - 1)
                    this.(varargin{i}) = varargin{i+1};
                end
                this = build(this);
            end
        end
        
        %% Parameterization Update Methods
        function this = build(this)            
            curves = makeNonintersecting([this.InputGeometry.Curves]);
            [X,Y]  = curves.getMidpoint;
            
         	nPlanes = numel(this.InputGeometry);
            nCurves = numel(curves);
            
            inPlane    = false(nCurves,nPlanes);
            onBoundary = false(nCurves,nPlanes);
            normalX    = zeros(nCurves,nPlanes);
            normalY    = zeros(nCurves,nPlanes);
            for i = 1:nPlanes
                [inPlane(:,i),...
                    onBoundary(:,i),normal] = this.InputGeometry(i).inOn(X,Y);
                normalX(:,i) = normal(:,1);
                normalY(:,i) = normal(:,2);
            end
            notInPlane = ~all(inPlane|onBoundary,2);
            
            innerProduct = bsxfun(@times,normalX,...
                                         reshape(normalX,nCurves,1,nPlanes))...
                          +bsxfun(@times,normalY,...
                                         reshape(normalY,nCurves,1,nPlanes));
            innerProduct = min(min(innerProduct,[],3),[],2);
            
            internalInterface = any(onBoundary,2) & (innerProduct < -sqrt(eps));
            
            midpointDistance  = sqrt( bsxfun(@minus,X,X.').^2 ...
                                       +bsxfun(@minus,Y,Y.').^2);
            areSameMidpoint   =   midpointDistance ...
                                < max(max(midpointDistance))*sqrt(eps);
                            
            duplicateInterfaces = any(tril(areSameMidpoint,-1),1).';
            
            remove      =  notInPlane | duplicateInterfaces | internalInterface;
            keep        = ~remove;
            
            curves   = curves(keep);
            rotation = this.vRotation;
            if rotation ~= 0
                curves = rotate(curves,'Rotation',rotation,...
                                       'Position',this.vPosition);
            end
            this.Curves = curves;
            this        = sortCurves(this);
%             [c1,c2,In1,In2,On1,On2,IP] = splitCurves(this);
%          	
%             In1 = In1|On1;
%             In2 = In2|On2;
%             
%             [n,~] = size(IP);
%             for i=1:n
%                 if IP(i,3)>0
%                     In1(IP(i,1))=true;
%                     In2(IP(i,2))=false;
%                 elseif IP(i,3)<0
%                     In1(IP(i,1))=false;
%                     In2(IP(i,2))=false;  
%                 end
%                 %In1(IP(i,1))=true;
%               	%In2(IP(i,2))=false;
%             end
%             
%             c = [c1(In1);c2(In2)];
%             
%             rotation    = this.vRotation;
%             if rotation ~= 0
%                 c = rotate(c,'Position',this.vPosition,...
%                             'Rotation',rotation);
%             end
%             
%             this.Curves 	= c;
%             this = sortCurves(this);
        end
        
        %% Compuational Geometry Methods
        
%         function [In,On,N] = inOn(this,xIn,yIn)
%             [In1,On1,N1] = this.Geometry1.inOn(xIn,yIn);
%             [In2,On2,N2] = this.Geometry2.inOn(xIn,yIn);
%             
%             In = In1&In2;
%             On = (On1&In2)|(On2&In1)|(On1&On2);
%             
%             assert(~any(In&On),'Intersection2D:inOn',...
%                   'Point is simultaneously in the region and on its boundary');
% 
%            	N = N1+N2;
%             N(~On,:) = 0;
%             N = N./repmat(sqrt(sum(N.^2,2)),1,2);
%             N(isnan(N))=0;
%         end
    end
end