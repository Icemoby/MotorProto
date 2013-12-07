conductorDiameter   = 1.0e-3;
insulationThickness = 0.1e-3;
nTurns              = 2;

%% Calculate x-coordinates which bounds the slot shape
outlineCurves       = conductorOutline.Curves;
xMin                = min([outlineCurves.vX0])*0.9;
xMax                = max([outlineCurves.vX0])*1.1;

%% Draw a line through the center of the slot
centerLine          = Geometry1D.draw('Polyline','Points',[xMin, 0; xMax, 0]);

%% Find where the center-line intersects the slot
[~,~,s]             = intersection([outlineCurves, centerLine]);
s                   = s(end);
x                   = centerLine.x(s{1});
xMin                = min(x);
xMax                = max(x);

%% Calculate maximum number of horizontal conductor layers on the center line
nHorzLayers         = floor((xMax - xMin) / (conductorDiameter + 2 * insulationThickness));

%% Determine centers of these conductors
xOffset             = conductorDiameter / 2 + (xMax - xMin - conductorDiameter * nHorzLayers) / (nHorzLayers + 1);
xCenter             = linspace(xMin + xOffset, xMax - xOffset, nHorzLayers);

%% Calculate y-coordinates which bound the slot shape
yMin                = min([outlineCurves.vY0])*1.1;
yMax                = max([outlineCurves.vY0])*1.1;

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
        
        in             = conductorOutline.inOn(xTest, yTest);
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
conductors  = Sector.empty(0,nConductors);
for i = 1:nConductors
    conductors(i) = Geometry2D.draw('Sector', 'Radius', [0 conductorDiameter / 2], 'Position', [xCenter(i),yCenter(i)], 'Angle', 2*pi);
end

% for i = 1:(nConductors / nTurns)
%     figure(2);clf;hold on;
%     outlineCurves.plot;
%     centerLine.plot;
%     verticalLines.plot;
%     scatter(xCenter,yCenter);
% 	conductors(seriesGroups(:,i)).wireframe;
%     pause;
% end
% 
% figure(2);clf;hold on;
% outlineCurves.plot;
% centerLine.plot;
% verticalLines.plot;
% scatter(xCenter,yCenter);
% for i = 1:nTurns
%     conductors(turnGroups(:,i)).wireframe;
%     pause;
% end