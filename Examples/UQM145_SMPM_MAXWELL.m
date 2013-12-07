%% IPM_Machine_Tutorial.m
clear all;
close all;

tic

%% Initialize the toolbox
simulation = MotorProto('IPM Machine Tutorial');
%warning on 'MotorProto:Verbose';

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('IPM Rotor','SynchronousRotor');
stator = model.newAssembly('IPM Stator','Stator');

%% Define General Machine Parameters
nPoles            = 18;
nTeethPerPhase    = 2;
nTeeth            = 3 * nTeethPerPhase * nPoles;
len               = 0.14;
statorOuterRadius = 0.125;
statorInnerRadius = 0.0982;
rotorOuterRadius  = 0.0972;
rotorInnerRadius  = 0.0834;
w_r               = 70;

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = len;
stator.Poles               = nPoles;
stator.Teeth               = 3 * nPoles * nTeethPerPhase;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = Steel1010;
stator.SourceType          = 'VoltageSource';
stator.ConnectionType      = 'Wye';
stator.ConductorDynamics   = 'Dynamic';
stator.Slot.Turns          = 3;

%Stranded-Style Conductors
% stator.Slot.ConductorType                 = 'Circular';
% stator.Slot.Conductor.ConductorDiameter   = 1.00e-3;
% stator.Slot.Conductor.InsulationThickness = 0.1e-3;

%Bus-Bar Style Conductors
stator.Slot.ConductorType           = 'Homogenized';
stator.Slot.Conductor.PackingFactor = 1;

%% Define slot geometry
slotWidth   = 0.5;
slotLength  = 0.7;
notchWidth  = 2 * 1 / 2 / nTeethPerPhase / 3 * 1.05;
notchLength = 1/29;

[slotOutline, ~] = slotTemplate(nTeeth, statorInnerRadius, statorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 1, 0);

conductorLength   = 13.5e-3;
conductorWidth    = 2.35e-3;
conductorPosition = [statorInnerRadius + conductorLength / 2 + 3.95e-3, 0];

conductorOutline = Geometry2D.draw('Rect', 'Width', conductorWidth, 'Length', conductorLength', 'Base', 'Center', 'Position', conductorPosition, 'PlotStyle', {'y'});

stator.Slot.Shape        = conductorOutline;
stator.ConductorMaterial = CopperExampleMaterial;

slotOutline = slotOutline - conductorOutline;
stator.addRegion('slot', slotOutline, Air, 'Static', -1); 

%% Set Rotor Parameters
rotor.Poles               = nPoles;
rotor.Length              = len;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = Steel1010;
rotor.OperatingMode       = 'synchronous'; %'locked'
rotor.InitialAngle        = 0;

%% Create Rotor Permanent Magnet
pmWidth    = 4e-3;
pmEmbrace  = (1-1/13)*1.02;
pmLength   = 2 * (rotorOuterRadius - pmWidth) * tan(2 * pi / nPoles / 2) * pmEmbrace;
pmPosition = [statorOuterRadius / 2 + rotorOuterRadius - pmWidth / 2, 0];

pmBody = Geometry2D.draw('Rect', 'Width', pmLength, 'Length', pmWidth + statorOuterRadius, 'Base', 'Center', 'Position', pmPosition, 'PlotStyle', {'m'});
pmTrim = Geometry2D.draw('Sector', 'Radius', [rotorInnerRadius, rotorOuterRadius], 'Angle', 2 * pi / nPoles, 'Rotation', - pi / nPoles);
permanentMagnet = pmBody * pmTrim;

rotor.addRegion('pm', permanentMagnet,  NdFe35, 'Dynamic', 0);

%% Trim Iron Between Magnets
trim1 = [rotorOuterRadius - pmWidth, pmLength / 2];
trim2 = [sqrt(statorOuterRadius^2 - pmLength^2 / 4), pmLength /2];
trim3 = statorOuterRadius * [cos(pi / nPoles), sin(pi / nPoles)];

poleM  = tan( pi / nPoles);
trim4M = - trim1(1) / trim1(2);
trim4b = trim1(2) - trim4M * trim1(1);

trim4 = trim4b / (poleM - trim4M) * [1, poleM];

trimPoints = [trim1;trim2;trim3;trim4];

trimUHP = Geometry2D.draw('Polygon2D', 'Points', trimPoints, 'PlotStyle', {'w'});
trimUHP = trimUHP * pmTrim;

trimPoints(:,2) = -trimPoints(:,2);
trimPoints      = flipud(trimPoints);

trimLHP = Geometry2D.draw('Polygon2D', 'Points', trimPoints, 'PlotStyle', {'w'});
trimLHP = trimLHP * pmTrim;

rotor.addRegion('trimUHP', trimUHP, Air, 'Static', 'None');
rotor.addRegion('trimLHP', trimLHP, Air, 'Static', 'None');

%% Set mesh parameters
mesh                       = simulation.Mesh;
mesh(1).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 160;
mesh(1).UseUniformGrid     = true;
mesh(2).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 160;
mesh(2).UseUniformGrid     = true;

%% Set Excitation
stator.Sources.ElectricalFrequency = w_r * nPoles / 2;

%Voltage Source
% stator.SourceType = 'VoltageSource';
% stator.Sources.HarmonicNumbers    = 1;
% stator.Sources.HarmonicAmplitudes = 708;
% stator.Sources.HarmonicPhases     = -2*pi*101/360;

%Current Source
stator.SourceType = 'CurrentSource';
stator.Sources.HarmonicNumbers    = 1;
stator.Sources.HarmonicAmplitudes = 150 / sqrt(3);
stator.Sources.HarmonicPhases     = pi * (1 + 1/6);

%% Configure algorithm
nSlotHarmonics = 2;
nTimePoints    = 2 * (2 * nTeethPerPhase * 3 + 1) * nSlotHarmonics + 1;
% simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'Verbose', true);
simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true, 'ShootingTolerance', 1e-3);
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true, 'TransientSolver', true, 'ShootingTolerance', 1e-3);
% simulation.configureAlgorithm('HarmonicBalance', 'TimePoints', nTimePoints, 'Verbose',true, 'NewtonTolerance', 1e-3);
% simulation.configureAlgorithm('IterativeHarmonicBalance', 'TimePoints', nTimePoints);
% simulation.configureAlgorithm('ParallelIterativeHarmonicBalance', 'TimePoints', nTimePoints); 
% % 
% simulation.configureAlgorithm('ParallelStatic', 'TimePoints', nTimePoints, 'Verbose', true);
solution = simulation.run;

% shutdownWorkers;
%% Plotting
% solution.plot('A','Time',1);
% solution.plot('B','Time',1);
% solution.plot('A','Harmonic',[0, 1]);
% % solution.plot('B','Harmonic',[0, 1]);
% solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% % solution.plot('J','Harmonic',1);
% % solution.plot('J','Time',1);
% % 
% % solution.plot('Flux Linkage','Time');
% % solution.plot('Flux Linkage','Harmonic');
solution.plot('Torque','Time');
% solution.plot('Torque','Harmonic');
% solution.plot('Voltage','Time');
% solution.plot('Voltage','Harmonic');
% solution.plot('Current','Time');
% solution.plot('Current','Harmonic');
% % 
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true);
% solution = simulation.run;
% solution.plot('Flux Linkage','Time');
% solution.plot('Flux Linkage','Harmonic');
% % solution.plot('Torque','Time');
% solution.plot('Voltage','Time');
% solution.plot('Voltage','Harmonic');
% % solution.plot('Current','Time');