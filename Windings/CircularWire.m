classdef CircularWire < Wire
    properties
        ConductorDiameter   = 0;
        InsulationThickness = 0;
    end
    
    methods
        %% Build
        function [conductors, nonConductors, connectionMatrix] = build(this, slotShape, conductorBoundaries, nTurns, nLayers, conductorDynamics, windingType, label)
            %% Get parameters
            conductorDiameter   = this.ConductorDiameter;
            insulationThickness = this.InsulationThickness;
            
            %% Calculate x-coordinates which bounds the slot shape
          	outlineCurves = slotShape.Curves;
            xMin          = min([outlineCurves.vX0])*0.9;
            xMax          = max([outlineCurves.vX0])*1.1;
            
            if windingType == WindingTypes.Concentrated
                yMin          = min([outlineCurves.vY0])*1.1;
                halfSlot      = Geometry2D.draw('Polygon2D','Points',[xMin,yMin;xMax,yMin;xMax,0;xMin,0]);
                halfSlot      = halfSlot * slotShape;
                outlineCurves = halfSlot.Curves;
          	end
            
            %% Draw a line through the center of the slot
            centerLine          = Geometry1D.draw('Polyline','Points',[xMin, 0; xMax, 0]);
            
            %% Find where the center-line intersects the slot
            [~,~,s]             = intersection([outlineCurves, centerLine]);
            s                   = s(end);
            x                   = centerLine.x(s{1});
            xMin                = min(x);
            xMax                = max(x);
            
            %% Restrict center-line to line between the conductorBoundaries
            xMin                = max(xMin, conductorBoundaries(1));
            xMax                = min(xMax, conductorBoundaries(2));
            
            %% Calculate maximum number of horizontal conductor layers on the center line
            nHorzLayers         = floor((xMax - xMin) / (conductorDiameter + 2 * insulationThickness));
            
            %% Determine centers of these conductors
            xOffset             = conductorDiameter / 2 + (xMax - xMin - conductorDiameter * nHorzLayers) / (nHorzLayers + 1);
            xCenter             = linspace(xMin + xOffset, xMax - xOffset, nHorzLayers);
            
            %% Calculate y-coordinates which bound the slot shape
           	yMin = min([outlineCurves.vY0]);
            yMax = max([outlineCurves.vY0]);
            
            dy   = (yMax - yMin)*0.1;
            
            yMin = yMin-dy;
            yMax = yMax+dy;
            
            %% Draw a set of verticle lines through the conductors on the center line
            verticalLines       = Polyline.empty(0,nHorzLayers);
            
            for i = 1:nHorzLayers
                verticalLines(i) = Geometry1D.draw('Polyline','Points',[xCenter(i), yMin; xCenter(i), yMax]);
            end
            
            %% Calculate intersection of the vertical lines with the slot and center-line
            [~,~,s] = intersection([outlineCurves, verticalLines]);
            s       = s(end-nHorzLayers+1:end);
            
            %% Determine number of vertical layers
            nVertLayers = zeros(1,nHorzLayers);
            yCenter     = cell(1,nHorzLayers);
            for i = 1:nHorzLayers
                %% First come up with a rough estimate of the number of vertical layers
                y              = verticalLines(i).y(s{i});
                yMin           = min(y);
                yMax           = max(y);
                nVertLayers(i) = floor((yMax - yMin) / (conductorDiameter + 2 * insulationThickness));
                
                %% Determine a set of test points
                in             = false;
                while ~all(in)
                    yOffset        = conductorDiameter / 2 + (yMax - yMin - conductorDiameter * nVertLayers(i)) / (nVertLayers(i) + 1);
                    yCenter{i}     = linspace(yMin + yOffset, yMax - yOffset, nVertLayers(i));
                    
                    rTest          = conductorDiameter / 2 + insulationThickness;
                    aTest          = linspace(0,2*pi*7/8,8).';
                    
                    yTest          = bsxfun(@plus, yCenter{i}, rTest * sin(aTest));
                    yTest          = reshape(yTest,[],1);
                    
                    xTest          = bsxfun(@plus, yCenter{i} * 0 + xCenter(i), rTest * cos(aTest));
                    xTest          = reshape(xTest,[],1);
                    
                    in             = slotShape.inOn(xTest, yTest);
                    if ~all(in)
                        nVertLayers(i) = nVertLayers(i) - 1;
                    end
                end
            end
            
            %% Ensure the number of conductors is divisible by the number of turns
            nConductors  = sum(nVertLayers);
            nRemoveTotal = mod(nConductors, nTurns);
            i = 1;
            while nRemoveTotal > 0
                if nVertLayers(i) < nRemoveTotal
                    nRemoveTotal   = nRemoveTotal - nVertLayers(i);
                    nVertLayers(i) = 0;
                    i = i + 1;
                else
                    nVertLayers(i) = nVertLayers(i) - nRemoveTotal;
                    nRemoveTotal   = 0;
                end
            end
            
            %% Adjust vertical layers
            xCenter     = num2cell(xCenter);
            for i = 1:nHorzLayers
                if nVertLayers(i) > 0
                    y              = verticalLines(i).y(s{i});
                    yMin           = min(y);
                    yMax           = max(y);
                    yOffset        = conductorDiameter / 2 + (yMax - yMin - conductorDiameter * nVertLayers(i)) / (nVertLayers(i) + 1);
                    yCenter{i}     = linspace(yMin + yOffset, yMax - yOffset, nVertLayers(i));
                    xCenter{i}     = ones(1,nVertLayers(i)) * xCenter{i};
                else
                    yCenter{i}     = [];
                    xCenter{i}     = [];
                end
            end
            
            %% Determine turn groups
            nConductors = sum(nVertLayers);
            yCenter     = cell2mat(yCenter);
            xCenter     = cell2mat(xCenter);
            
            assert(nConductors > 0,'MotorProto:CircularWire','Not able to place conductors in slot. Check conductor size');
            
            xcMin       = min(xCenter);
            ycMin       = min(yCenter);
            
            dX          = xCenter - xcMin;
            dY          = yCenter;
            
            dR          = hypot(dX,dY);
            dA          = atan2(dY,dX);
            dA          = round(dA / sqrt(eps)) * sqrt(eps);
            
            [~,I]       = sortrows([dR;dA].',[1 2]);
            
            turnGroups  = reshape(I,[],nTurns);
            
            %% Determine series groups
            seriesGroups      = turnGroups.';
            [~,J]             = sort(xCenter(seriesGroups(1,:)) + sqrt(eps) * abs(yCenter(seriesGroups(1,:))), 'descend');
            seriesGroups(1,:) = seriesGroups(1,J);
            for i = 2:nTurns
                [~,J] = sort(xCenter(seriesGroups(i-1,:)) + sqrt(eps) * abs(yCenter(seriesGroups(i-1,:))), 'descend');
                I     = 1:numel(J);
                
                meritMatrix = abs(bsxfun(@minus, xCenter(seriesGroups(i-1,:)), xCenter(turnGroups(:,i)).')) - sqrt(eps) * abs(bsxfun(@minus, yCenter(seriesGroups(i-1,:)), yCenter(turnGroups(:,i)).'));
                for j = 1:numel(J)
                    [~,k]                = max(meritMatrix(I, J(j)));
                    seriesGroups(i,J(j)) = turnGroups(I(k),i);
                    I(k)                 = [];
                end
            end
            
            %% Draw conductors
            switch windingType
                case WindingTypes.Distributed
                    conductorGeometry = Sector.empty(0,nConductors);
                case WindingTypes.Concentrated
                    conductorGeometry = Sector.empty(0,2*nConductors);
                otherwise
                    error('');
            end
            
            for i = 1:nConductors
                conductorGeometry(i) = Geometry2D.draw('Sector', 'Radius', [0 conductorDiameter / 2], 'Position', [xCenter(i),yCenter(i)], 'Angle', 2*pi,'PlotStyle',{'y'});
            end
            
            switch windingType
                case WindingTypes.Concentrated
                    for i = 1:nConductors
                        conductorGeometry(i+nConductors) = Geometry2D.draw('Sector', 'Radius', [0 conductorDiameter / 2], 'Position', [xCenter(i),-yCenter(i)], 'Angle', 2*pi,'PlotStyle',{'y'});
                    end
                    nConductors = nConductors * 2;
                case WindingTypes.Distributed
                    %NOOP
                otherwise
                    error('');
            end
            
            %% Assign Regions
            nonConductors = Region2D('Geometry', slotShape - conductorGeometry, 'Material', this.InsulatorMaterial, 'Dynamics', DynamicsTypes.Static,...
                                     'Name', [label , ' NC']);
            
            conductors = Region2D.empty(0,nConductors);
            for i = 1:nConductors
                conductors(i) = Region2D('Geometry', conductorGeometry(i), 'Material', this.ConductorMaterial, 'Dynamics', conductorDynamics,...
                                         'Name', [label , ' C', num2str(i)]);
            end
            
            connectionMatrix = seriesGroups.';
            connectionMatrix = sortrows(connectionMatrix,1:nTurns);
            
            switch windingType
                case WindingTypes.Distributed
                    %NOOP
                case WindingTypes.Concentrated
                    connectionMatrix = cat(3,connectionMatrix,connectionMatrix+nConductors/2);
                otherwise
                    error('');
            end
        end
    end
end