classdef HomogenizedConductor < Wire
    %HomogenizedConductor.m A class representing a bulk averaged conductor model
    %   HomogenizedConductor objects
    %
    % HomogenizedConductor properties:
    %   PackingFactor       - Ratio of conducting material to total material
    %   HomogenizedMaterial - Output of the material homogenization process
    %
  	% HomogenizedConductor methods:
    %   build - Divides an input region into a number of conductors
    %
    % HomogenizedConductor inherits properties from Wire.
    %
    % See the help for Wire for more information.
    %
    % See also MotorProto, Model, Assembly, Wire
    
    %   Copyright 2012 Jason Pries
    %   $Revision 0.0.0.1$
    
%{
properties:
 	%PackingFactor - Ratio of conducting material to total material
    %   
    %
    %
    % See also HomogenizedConductor
    PackingFactor;
    
 	%HomogenizedMaterial - Output of the material homogenization process
    %
    % See also HomogenizedConductor, HomogenizedMaterialProperty
    HomogenizedMaterial;
%}
    
    properties
        PackingFactor = HomogenizedConductor.setProperty(1);
    end
    
    properties (Dependent)
        HomogenizedMaterial
    end
    
    methods
        %% Setters
        function this = set.PackingFactor(this, value)
            this.PackingFactor = this.setProperty(value);
        end
        
        %% Getters
        function value = get.PackingFactor(this)
            value = this.PackingFactor.Value;
        end
        
        function value = get.HomogenizedMaterial(this)
            m     = [this.ConductorMaterial, this.InsulatorMaterial];
            p     = [this.PackingFactor, 1 - this.PackingFactor];
            value = HeterogeneousMaterial('BaseProperties', m, 'Percentage', p);
        end
        
        %% Build
        function [conductors, nonConductors, connectionMatrix] = build(this, slotShape, conductorBoundaries, nTurns, nLayers, conductorDynamics, windingType, label)
           	 %% Turn Division by Sequential-Quadratic-Interpolation
            
            %% Get Bounding Box of slotShape
            curves  = slotShape.Curves;
            nCurves = numel(curves);
            xBounds = [inf 0];
            yBounds = [inf -inf];
            for i = 1:nCurves;
                xBounds(1)  = min(xBounds(1), curves(i).bbCenter(1) - curves(i).bbRadius);
                xBounds(2)  = max(xBounds(2), curves(i).bbCenter(1) + curves(i).bbRadius);
                
                yBounds(1)  = min(yBounds(1), curves(i).bbCenter(2) - curves(i).bbRadius);
                yBounds(2)  = max(yBounds(2), curves(i).bbCenter(2) + curves(i).bbRadius);
            end
            
            xMin = 0;
            yMin = yBounds(1);
            yMax = yBounds(2);
            
            %% Divide into Conducting/Nonconducting Regions
            if ~isinf(conductorBoundaries(1)) && ~isinf(conductorBoundaries(2))
                %% Two Boundaries
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [conductorBoundaries(1), yMin;
                                                                                  conductorBoundaries(2), yMin;
                                                                                  conductorBoundaries(2), yMax;
                                                                                  conductorBoundaries(1), yMax]);
                conductingRegion        = slotShape * testShape;
                
                nonConductors           = Region2D.empty(0,2);
                
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [conductorBoundaries(1), yMin;
                                                                                  conductorBoundaries(1), yMax;
                                                                                  xBounds(1)            , yMax;
                                                                                  xBounds(1)            , yMin]);
                                                                              
                nonConductors(1)        = Region2D('Geometry', slotShape * testShape, 'Material', this.InsulatorMaterial, 'Dynamics', DynamicsTypes.Static,...
                                                   'Name', [label, ' NC1']);
                
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [conductorBoundaries(2), yMin;
                                                                                  xBounds(2)            , yMin;
                                                                                  xBounds(2)            , yMax;
                                                                                  conductorBoundaries(2), yMax]);
                                                                              
                nonConductors(2)        = Region2D('Geometry', slotShape * testShape, 'Material', this.InsulatorMaterial, 'Dynamics', DynamicsTypes.Static,...
                                                   'Name', [label, ' NC2']);
                
                xBounds(1)              = conductorBoundaries(1) * (1 - sqrt(eps));
                xBounds(2)              = conductorBoundaries(2) * (1 + sqrt(eps));
            elseif ~isinf(conductorBoundaries(1))
                %% Inner Boundary Only
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [conductorBoundaries(1), yMin;
                                                                                  xBounds(2), yMin;
                                                                                  xBounds(2), yMax;
                                                                                  conductorBoundaries(1), yMax]);
                conductingRegion        = slotShape * testShape;
                
                testShape               = Geometry2D.draw('Polygon2D','Points', [conductorBoundaries(1), yMin;
                                                                                 conductorBoundaries(1), yMax;
                                                                                 xBounds(1)            , yMax;
                                                                                 xBounds(1)            , yMin]);
                nonConductors           = Region2D('Geometry', slotShape * testShape, 'Material', this.InsulatorMaterial, 'Dynamics', DynamicsTypes.Static,...
                                                   'Name', [label, ' NC']);
                
                xBounds(1)              = conductorBoundaries(1) * (1 - sqrt(eps));
            elseif ~isinf(conductorBoundaries(2))
                %% Outer Boundary Only
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [xBounds(1), yMin;
                                                                                  conductorBoundaries(2), yMin;
                                                                                  conductorBoundaries(2), yMax;
                                                                                  xBounds(1), yMax]);
                conductingRegion        = slotShape * testShape;
                
                testShape               = Geometry2D.draw('Polygon2D', 'Points', [conductorBoundaries(2), yMin;
                                                                                  xBounds(2)            , yMin;
                                                                                  xBounds(2)            , yMax;
                                                                                  conductorBoundaries(2), yMax]);

                nonConductors           = Region2D('Geometry', slotShape * testShape, 'Material', this.InsulatorMaterial, 'Dynamics', DynamicsTypes.Static,...
                                                   'Name', [label, ' NC']);
                
                xBounds(2)              = conductorBoundaries(2) * (1 + sqrt(eps));
            else
                %% No Boundaries
                conductingRegion = slotShape;
                nonConductors    = Region2D.empty(0,1);
            end
            
            %% Get Good Initial Polynomial Area Estimate
           	A  = area(conductingRegion);
            dA = A / nTurns;
            x  = [(xBounds(1)*3+xBounds(2))/4 (xBounds(1)+xBounds(2))/2 (xBounds(1)+xBounds(2)*3)/4];
            a  = zeros(1,3);
            for i = 1:3
                testRegion = Geometry2D.draw('Polygon2D','Points', [xMin, yMin;x(i), yMin; x(i), yMax; xMin, yMax]);
                testTurn   = slotShape * testRegion;
                a(i)       = area(testTurn);
            end
            
            %% Divide Conducting Region in Massive Turns
            conductorGeometry = Polygon2D.empty(0, nTurns);
            for i = 1:nTurns
                aNew = 0;
                while abs(aNew - dA) > dA * sqrt(eps)
                    p      = polyfit(x,a,2);
                    p(end) = p(end) - i*dA;
                    
                    if i == nTurns
                        xMax = max(xBounds) * 2;
                    else
                        xMax = real(roots(p));

                        isInRange = (xMax > xBounds(1)) & (xMax < xBounds(2));
                        if any(isInRange)
                            xMax = max(xMax(isInRange));
                        else
                            xMax = xBounds(2);
                        end
                    end
                    
                    testRegion           = Geometry2D.draw('Polygon2D', 'Points', [xMin, yMin; xMax, yMin; xMax, yMax; xMin, yMax]);
                    conductorGeometry(i) = conductingRegion * testRegion;
                    aNew                 = area(conductorGeometry(i));
                    
                    if abs(aNew - dA) > dA * sqrt(eps)
                        x = [x(2:3) xMax];
                        a = [a(2:3) (i-1)*dA + aNew];
                    end
                end
                xMin = xMax;
            end
            
            if (nLayers == 2) && (windingType == WindingTypes.Concentrated)
                conductors = Region2D.empty(0, 2 * nTurns);
                testRegion = Geometry2D.draw('Polygon2D','Points',[xBounds(1), 0; xBounds(2), 0; xBounds(2), yBounds(2); xBounds(1), yBounds(2)]);
                
                for i = 1:nTurns
                    conductors(i) = Region2D('Geometry', conductorGeometry(i) * testRegion, 'Material', this.HomogenizedMaterial, 'Dynamics', conductorDynamics,'Name', [label, ' C', num2str(i)]);
                end
                
                for i = (nTurns+1):(2*nTurns)
                    conductors(i) = Region2D('Geometry', conductorGeometry(i-nTurns) - testRegion, 'Material', this.HomogenizedMaterial, 'Dynamics', conductorDynamics,'Name', [label, ' C', num2str(i)]);
                end
                
                connectionMatrix        = zeros(nTurns,1,2);
                connectionMatrix(:,1,1) = (1:nTurns)';
                connectionMatrix(:,1,2) = ((nTurns+1):(2*nTurns))';
            elseif nLayers == 1
                conductors = Region2D.empty(0, nTurns);
                for i = 1:nTurns
                    conductors(i) = Region2D('Geometry', conductorGeometry(i), 'Material', this.HomogenizedMaterial, 'Dynamics', conductorDynamics,'Name', [label, ' C', num2str(i)]);
                end
                
             	connectionMatrix = (1:nTurns)';
            else
                error('Number of winding Layers must be equal to 1 or 2');
            end
        end
    end
end