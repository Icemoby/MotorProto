%% IPM_Machine_Tutorial.m
clear all;
close all;
tic

%% Initialize the toolbox
simulation = MotorProto('IPM Machine Tutorial');
%warning on 'MotorProto:Verbose';

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('IPM Rotor','InductionRotor');
stator = model.newAssembly('IPM Stator','Stator');

%% Define General Machine Parameters
nPoles            = 8;
nTeethPerPhase    = 3;
nTeeth            = 3 * nTeethPerPhase * nPoles;
len               = 0.18;
statorOuterRadius = 0.12;
statorInnerRadius = 0.08513;
rotorOuterRadius  = 0.085;
rotorInnerRadius  = 0.05;
w_r               = 70;

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = len;
stator.Poles               = nPoles;
stator.Teeth               = 3 * nPoles * nTeethPerPhase;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = Arnon5;
stator.SourceType          = 'CurrentSource';
stator.ConnectionType      = 'Wye';
stator.ConductorDynamics   = 'Dynamic';
stator.Slot.Turns          = 2;

% % Stranded-Style Conductors
% stator.Slot.ConductorType                 = 'Circular';
% stator.Slot.Conductor.ConductorDiameter   = 1.75e-3;
% stator.Slot.Conductor.InsulationThickness = 0.1e-3;

% Bus-Bar Style Conductors
stator.Slot.ConductorType           = 'Homogenized';
stator.Slot.Conductor.PackingFactor = 0.5;

%% Define slot geometry
slotWidth   = 0.5;
slotLength  = 0.5;
notchWidth  = 1 / 2 / nTeethPerPhase / 3;
notchLength = 0.01;

[conductorOutline, innerSlot] = slotTemplate(nTeeth, statorInnerRadius, statorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 0, 'Auto');
stator.Slot.Shape             = conductorOutline;
stator.ConductorMaterial      = CopperExampleMaterial;
                
stator.addRegion('innerSlot', innerSlot, Air, 'Static', -1); 

%% Set Rotor Parameters
rotor.Poles               = nPoles;
rotor.Teeth               = (3 * nTeethPerPhase + 1) * nPoles;
rotor.Length              = len;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.Slip                = 0.01;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = Arnon5;
rotor.InitialAngle        = 0;

%% Draw Key-Style Rotor Slot
conductorRadius = 0.25 * (2 * rotorOuterRadius * sin(pi / (3 * nTeethPerPhase + 1) / nPoles));
conductorInset  = 2    * (statorInnerRadius - rotorOuterRadius);
notchAngle      = 2    * pi / ((3 * nTeethPerPhase + 1) * nPoles) * 0.1;

slotNotch = Geometry2D.draw('Sector','Radius',[rotorOuterRadius - conductorRadius / 2 - conductorInset, rotorOuterRadius],...
                                     'Angle' , notchAngle, 'Rotation', -notchAngle / 2);

slotShape = Geometry2D.draw('Sector','Radius',[0, conductorRadius],'Angle',2*pi,'Position',[rotorOuterRadius - conductorRadius - conductorInset, 0]);

slotNotch = slotNotch - slotShape;

rotor.addRegion('slotShape', slotShape, CopperExampleMaterial, 'Floating', 'None');
rotor.addRegion('slotNotch', slotNotch, Air, 'Static', 'None');

%% Set mesh parameters
mesh                       = simulation.Mesh;
mesh(1).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 40;
mesh(2).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 40;

%% Set Excitation
stator.Sources.ElectricalFrequency = w_r * nPoles / 2;

%Voltage Source
% stator.Sources.HarmonicNumbers    = 1:2:17;
% stator.Sources.HarmonicAmplitudes = [392.940960259865,38.7974461566293,21.1686231750374,18.5295847823860,6.54971559669156,2.95498716209424,8.02036987709044,4.85090773859384,6.58391266174923;];
% stator.Sources.HarmonicPhases     = [1.83559893815957,2.83723513902788,-2.53101267526780,-2.54878725386589,-3.09621299590694,1.63134692441761,0.313394958242182,-1.28085787831664,2.22475806066111;];

%Current Source
stator.Sources.HarmonicNumbers    = 1;
stator.Sources.HarmonicAmplitudes = 305 / sqrt(3) * 0.01;
stator.Sources.HarmonicPhases     = (2 * pi / 3) + 2 * pi / 3 * 1.1;

%% Configure algorithm
nSlotHarmonics = 1;
nTimePoints    = 2 * (2 * nTeethPerPhase * 3 * + 1) * nSlotHarmonics + 1;

% simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'Verbose', true);
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true);
simulation.configureAlgorithm('StroboscopicShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true);
% simulation.configureAlgorithm('HarmonicBalance', 'TimePoints', nTimePoints);
% error('Don't use this'); simulation.configureAlgorithm('IterativeHarmonicBalance', 'TimePoints', nTimePoints, 'ReportProgress', true); 

solution = simulation.run;

%% Plotting
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