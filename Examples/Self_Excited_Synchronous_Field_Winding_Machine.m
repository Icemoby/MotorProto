%% UQM145_SMPM_Tutorial.m
clear all;
close all;

tic

%% Initialize the toolbox
simulation = MotorProto('UQM145_SMPM');
%warning on 'MotorProto:Verbose';

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('SESFW Rotor','SelfExcitedSynchronousRotor');
stator = model.newAssembly('SESFW Stator','Stator');

%% Define General Machine Parameters
cx = 0.141643;
cy = 0.156408;
nPoles = 4;
nStatorTeethPerPhase = 3;
nStatorTeeth = 3 * nPoles * nStatorTeethPerPhase;
stackLength = 0.128;
statorOuterRadius = 0.193/2;
statorInnerRadius = 0.113514/2;
rotorOuterRadius = 0.112523/2;
rotorInnerRadius = 0.030/2;
w_r = 25;

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = stackLength;
stator.Poles               = nPoles;
stator.Teeth               = nStatorTeeth;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = Steel1010;
stator.SourceType          = 'CurrentSource';
stator.ConnectionType      = 'Wye';
stator.ConductorDynamics   = 'Dynamic';
stator.WindingType         = 'Distributed';
stator.Slot.Turns          = 14;

%% Stranded-Style Conductors
stator.Slot.ConductorType                 = 'Circular';
stator.Slot.Conductor.ConductorDiameter   = 0.001016*1.8;
stator.Slot.Conductor.InsulationThickness = 0.0001016*1.8;
% stator.Slot.ConductorType           = 'Homogenized';
% stator.Slot.Conductor.PackingFactor = 0.5;

%% Define stator slot geometry
backIronLength = 1-(0.153467-0.113514)/(0.193-0.113514);
slotLength = (0.016182+0.002+0.001)/((0.153467-0.113514)/2);
notchLength = 1 - slotLength;
slotLength = slotLength * (1 - backIronLength);
notchLength = notchLength * (1 - backIronLength);

slotWidth   = 0.5;
notchWidth  = 1.525/9.894;

[slotOutline, slotNotch] = slotTemplate(nStatorTeeth, statorInnerRadius, statorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 0.5, 0,'InnerSlotShape','rounded','OuterSlotShape','rounded');

stator.Slot.Shape        = slotOutline;
stator.ConductorMaterial = CopperExampleMaterial;

if ~isempty(slotNotch)
    stator.addRegion('slot', slotNotch, Air, 'Static', -1); 
end

%% Set Rotor Parameters
nRotorFieldSlots = 5 * 4;
nRotorTransformerSlots = 2 * 4;
nRotorTeeth = nRotorFieldSlots + nRotorTransformerSlots;

rotor.Poles               = nPoles;
rotor.Teeth               = nRotorTeeth;
rotor.FieldSlots          = nRotorFieldSlots;
rotor.TransformerSlots    = nRotorTransformerSlots;
rotor.Length              = stackLength;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = Steel1010;
rotor.OperatingMode       = 'locked';
rotor.InitialAngle        = (6/18-1+2/14)*pi/2;
rotor.ConductorDynamics   = 'Dynamic';
rotor.Slot.Turns          = 45;

%%Stranded-Style Conductors
rotor.Slot.ConductorType                 = 'Circular';
rotor.Slot.Conductor.ConductorDiameter   = 0.001016;
rotor.Slot.Conductor.InsulationThickness = 0.0001016;
% rotor.Slot.Conductor.ConductorDiameter   = 0.00082;
% rotor.Slot.Conductor.InsulationThickness = 0.000082;

%% Define rotor slot geometry
backIronPercentage = (80.385-30)/(112.521-30);
notchLength = (0.563 * 2.0) / ((112.521-80.385) / 2);
slotLength = 1 - notchLength;
notchLength = notchLength * (1 - backIronPercentage);
slotLength = slotLength * (1 - backIronPercentage);

slotWidth   = 7.187/(7.187+4.788);
notchWidth  = 2.805/(2.805+9.807);

[slotOutline, slotNotch] = slotTemplate(nRotorTeeth, rotorInnerRadius, rotorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 1, 0.033, 'InnerSlotShape','rounded','OuterSlotShape','rounded','AirgapLocation','outside');

rotor.Slot.Shape        = slotOutline;
rotor.ConductorMaterial = CopperExampleMaterial;

rotor.addRegion('slot', slotNotch, Air, 'Static', -1); 

%% Set mesh parameters
mesh                       = simulation.Mesh;
mesh(1).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 1;
mesh(2).MaximumElementSize = (statorOuterRadius - rotorInnerRadius) / 1;
mesh(2).UseUniformGrid = true;

%% Set Excitation
stator.Sources.ElectricalFrequency = w_r * nPoles / 2;
% rotor.Sources.ElectricalFrequency = w_r * nPoles / 2;

%% Current Source
stator.Sources.HarmonicNumbers    = 1;
stator.Sources.HarmonicAmplitudes = 0.8165 / sqrt(3);
stator.Sources.HarmonicPhases     = 0;

%% Configure algorithm
nSlotHarmonics = 1;
nTimePoints    = 2 * (2 * nStatorTeethPerPhase * 3 + 1) * nSlotHarmonics + 1;
% simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'Verbose', true);
simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true);
% simulation.configureAlgorithm('TPFEM', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true, 'MaxGMRESIterations', 10);
solution = simulation.run;

%% Plotting
solution.plot('A','Time',1);
solution.plot('B','Time',1);
% % % solution.plot('H','Time',1);
% % % solution.plot('M','Time',1);
% % % % solution.plot('A','Harmonic',[0, 1]);
solution.plot('B','Harmonic',[1]);
% solution.plot('H','Harmonic',[1]);
% % solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% solution.plot('J','Harmonic',1);
% % solution.plot('J','Time',1);
% % 
solution.plot('Flux Linkage','Time');
% solution.plot('Flux Linkage','Harmonic');
% % solution.plot('Torque','Time');
% % % solution.plot('Torque','Harmonic');
solution.plot('Voltage','Time');
solution.plot('Voltage','Harmonic');
solution.plot('Current','Time');
solution.plot('Current','Harmonic');