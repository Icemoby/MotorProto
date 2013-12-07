%% IPM_Machine_Tutorial.m
clear all;
close all;
tic
display('IPM_MACHINE_TEST_SCRIPT');

%% Initialize the toolbox
simulation = MotorProto('Interior Permanent Magnet Machine Demo');
%warning on 'MotorProto:Verbose'

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('IPM Rotor','SynchronousRotor');
stator = model.newAssembly('IPM Stator','Stator');

%% Define IPM Rotor Geometry and material properties
parameters = simulation.Parameters;
parameters.new('nPoles',            8);
parameters.new('nTeethPerPhase',	3);
parameters.new('len',               0.18);
parameters.new('statorOuterRadius', 0.12);
parameters.new('statorInnerRadius', 0.08513);
parameters.new('rotorOuterRadius',  0.085);
parameters.new('rotorInnerRadius',  0.05);
parameters.new('slotInsetPercent',  0.02);
parameters.new('slotWidthPercent',  0.5);
parameters.new('slotLengthPercent', 0.5);
parameters.new('toothGapPercent',   0.1);
parameters.new('w_r',               70);

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = len;
stator.Poles               = nPoles;
stator.Teeth               = 3 * nPoles * nTeethPerPhase;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = IronExampleMaterial;
stator.SourceType          = 'CurrentSource';
stator.ConnectionType      = 'Wye';

stator.ConductorDynamics                  = 'Dynamic';
stator.Slot.Turns                         = 2;
stator.Slot.ConductorType                 = 'Circular';
% stator.Slot.Conductor.ConductorDiameter   = 1e-3;
stator.Slot.Conductor.ConductorDiameter   = 3e-3;
stator.Slot.Conductor.InsulationThickness = 0.1e-3;

%% Define rough slot geometry
toothWidthAngle   = 2 * pi / 3 / nPoles / nTeethPerPhase;
slotInnerRadius   =   statorInnerRadius * (1 - slotInsetPercent)...
                    + statorOuterRadius * slotInsetPercent;
slotInnerPosition = [slotInnerRadius, 0];
slotInnerWidth    = slotInnerRadius * sin(toothWidthAngle * slotWidthPercent);
slotLength        = (statorOuterRadius - statorInnerRadius) * slotLengthPercent;
toothWidthAngle   = 0;

%% First create a polygon
slotPoints = [slotInnerPosition(1),...
              slotInnerPosition(2) + slotInnerWidth / 2 ;
              slotInnerPosition(1),...
              slotInnerPosition(2) - slotInnerWidth / 2 ];
slotPoints = [slotPoints;
              slotPoints(2,1) + slotLength * cos(toothWidthAngle / 2),...
              slotPoints(2,2) - slotLength * sin(toothWidthAngle / 2);
              slotPoints(1,1) + slotLength * cos(toothWidthAngle / 2),...
              slotPoints(1,2) + slotLength * sin(toothWidthAngle / 2)];
                   
conductorOutline = Geometry2D.draw('Polygon2D',...
                                    'Points',       slotPoints,...
                                    'PlotStyle',    {'y'});

%% Round outer slot edges
distanceSlotBack	= (slotPoints(4,2) - slotPoints(3,2)) / 2;
adjustmentFactor    = cos(toothWidthAngle / 2)...
                     +tan(toothWidthAngle / 2) * (1 + sin(toothWidthAngle/2));
roundingInnerRadius	= distanceSlotBack / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotBack^2);
roundingPosition	= [slotPoints(3,1) - roundingInnerRadius, 0];
roundingAngle       = pi + toothWidthAngle;
roundingRotation	= -(pi + toothWidthAngle) / 2;

outerRoundingSector = Geometry2D.draw('Sector',   ...
                                        'Radius',   [roundingInnerRadius,...
                                                        roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);
                                
conductorOutline	= conductorOutline - outerRoundingSector;

%% Round inner slot edges
distanceSlotFront	= slotInnerWidth / 2;
adjustmentFactor    =  cos(toothWidthAngle / 2)...
                      -tan(toothWidthAngle / 2) * (1 - sin(toothWidthAngle / 2));
roundingInnerRadius	= distanceSlotFront / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotFront^2);
roundingPosition	= [slotPoints(1,1) + roundingInnerRadius, 0];
roundingAngle       = pi - toothWidthAngle;
roundingRotation	= ( pi + toothWidthAngle ) / 2;

innerRoundingSector = Geometry2D.draw('Sector',   ...
                                      	'Radius',   [0, roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);

conductorOutline	= conductorOutline - innerRoundingSector;

innerSlot = Geometry2D.draw('Sector',...
                                'Radius',   [0,roundingInnerRadius],...
                                'Rotation', roundingRotation,...
                                'Angle',    roundingAngle,...
                                'Position', roundingPosition,...
                                'PlotStyle',{'w'});

%% Add gap between the teeth
gapInnerRadius = statorInnerRadius;
gapOuterRadius = statorInnerRadius + roundingInnerRadius;
gapWidthAngle  = 2 * pi / 3 / nPoles / nTeethPerPhase * toothGapPercent;

gap            = Geometry2D.draw('Sector',   ...
                                 	'Radius',	[gapInnerRadius,...
                                                	gapOuterRadius],...
                                  	'Angle',  	gapWidthAngle,...
                                	'Rotation',	-gapWidthAngle/2);
innerSlot      = innerSlot + gap;

stator.Slot.Shape        = conductorOutline;
stator.ConductorMaterial = CopperExampleMaterial;
                
stator.addRegion('innerSlot', innerSlot, Air, 'Static', -1); 

rotor.Poles               = nPoles;
rotor.Length              = len;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = IronExampleMaterial;
rotor.OperatingMode       = 'synchronous';
rotor.InitialAngle        = 0;

%% Create Rotor Permanent Magnet
parameters.new('pmAngle'         , 2 * pi / 3 / nPoles * 2);
% parameters.new('pmAngle'         , 2 * pi / nPoles * (2 * nTeethPerPhase * 3 - 1) / (2 * nTeethPerPhase * 3));
parameters.new('pmRadPercent'    , 0.9);
parameters.new('fbMinBridgeWidth', 0.0001);

dMax       = (rotorOuterRadius - fbMinBridgeWidth ) * cos(pmAngle / 2);
pmVolume   = 22 * (statorInnerRadius - rotorOuterRadius) * 2 * dMax * tan(pmAngle / 2);
dMin       = (rotorInnerRadius + fbMinBridgeWidth) * cos(pmAngle / 2);
d          = dMin - (dMin - dMax) * pmRadPercent;
pmLength   = 2 * d * tan(pmAngle / 2);
pmWidth    = pmVolume / pmLength;
pmPosition = [d - pmWidth / 2, 0];

permanentMagnet = Geometry2D.draw('Rect', 'Width', pmLength, 'Length', pmWidth, 'Base', 'Center', 'Position', pmPosition, 'PlotStyle',{'m'});
rotor.addRegion('pm', permanentMagnet,  PermanentMagnetExampleMaterial, 'Dynamic', 0);

%% Create Rotor Flux Barriers
rFB   = rotorOuterRadius - fbMinBridgeWidth;
pmRad = sqrt(d^2 + pmLength^2 / 4);

xFB0 = pmRad * cos(pmAngle / 2);
yFB0 = pmRad * sin(pmAngle / 2);

xFB1 = rFB * cos(pmAngle / 2);
yFB1 = rFB * sin(pmAngle / 2);

xFB2 = rFB * cos(pmAngle / 2 + pi * 1 / (3 * nTeethPerPhase * nPoles));
yFB2 = rFB * sin(pmAngle / 2 + pi * 1 / (3 * nTeethPerPhase * nPoles));

% xFB2 = rFB * cos(2 * pi / nPoles / 2);
% yFB2 = rFB * sin(2 * pi / nPoles / 2);

xFB4 = xFB0 - pmWidth;
yFB4 = yFB0;

rFB3 = sqrt(xFB4^2 + yFB4^2);
xFB3 = rFB3 * cos(pmAngle / 2 + pi * 1 / (3 * nTeethPerPhase * nPoles));
yFB3 = rFB3 * sin(pmAngle / 2 + pi * 1 / (3 * nTeethPerPhase * nPoles));
% xFB3 = rFB3 * cos(2 * pi / nPoles / 2);
% yFB3 = rFB3 * sin(2 * pi / nPoles / 2);


fbPoints = [xFB0 yFB0;
            xFB1 yFB1;
            xFB2 yFB2;
            xFB3 yFB3;
            xFB4 yFB4];
        
fb1 = Geometry2D.draw('Polygon2D', 'Points', fbPoints, 'PlotStyle', {'w'});
rotor.addRegion('fb1', fb1, Air, 'Static', -1);

fbPoints = [xFB4 -yFB4;
            xFB3 -yFB3;
            xFB2 -yFB2;
            xFB1 -yFB1;
            xFB0 -yFB0];
        
fb2 = Geometry2D.draw('Polygon2D', 'Points', fbPoints, 'PlotStyle', {'w'});
rotor.addRegion('fb2', fb2, Air, 'Static', -1);
        
%% Set mesh parameters
mesh                      = simulation.Mesh;
[mesh.MaximumElementSize] = deal((statorOuterRadius - rotorInnerRadius) / 40);

%% Set Excitation
% stator.SourceType         = 'voltage';
% stator.HarmonicNumbers    = 1:2:17;
% stator.HarmonicAmplitudes = [392.940960259865,38.7974461566293,21.1686231750374,18.5295847823860,6.54971559669156,2.95498716209424,8.02036987709044,4.85090773859384,6.58391266174923;];
% stator.HarmonicPhases     = [1.83559893815957,2.83723513902788,-2.53101267526780,-2.54878725386589,-3.09621299590694,1.63134692441761,0.313394958242182,-1.28085787831664,2.22475806066111;];

% stator.ConductorDynamics  = 'dynamic';
% stator.SourceType         = 'current';
% stator.HarmonicNumbers    = 1;
% stator.HarmonicAmplitudes = 290;
% stator.HarmonicPhases     = - pi / 2 + pi / 20;

stator.Sources.ElectricalFrequency = w_r * nPoles / 2;
% Wye Connected
stator.Sources.HarmonicAmplitudes  = 305 / sqrt(3);
stator.Sources.HarmonicPhases      = (2 * pi / 3) + 2 * pi / 3 * 1.1;
% % Delta Connected
% stator.Sources.HarmonicAmplitudes  = 305;
% stator.Sources.HarmonicPhases      = (2 * pi / 3) + 2 * pi / 3 * 1.1;

% stator.Sources.ElectricalFrequency = 280; 
% stator.Sources.HarmonicAmplitudes  = 393 / sqrt(3);
% stator.Sources.HarmonicPhases      = 1.84;

%% Configure algorithm
nTimePoints = 38*1+1;
simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'ReportProgress', true);
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true);
% simulation.configureAlgorithm('HarmonicBalance', 'TimePoints', nTimePoints);
% simulation.configureAlgorithm('ExplicitErrorCorrection', 'TimePoints', nTimePoints);

solution = simulation.run;
% 
% solution.plot('A','Time',1);
% solution.plot('B','Time',1);
% solution.plot('A','Harmonic',[0, 1]);
% solution.plot('B','Harmonic',[0, 1]);
% solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% solution.plot('J','Harmonic',1);
% solution.plot('J','Time',1);
% 
% solution.plot('Flux Linkage','Time');
% solution.plot('Flux Linkage','Harmonic');
% solution.plot('Torque','Time');
% solution.plot('Voltage','Time');
% solution.plot('Current','Time');