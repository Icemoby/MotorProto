%% UQM145_SMPM_Tutorial.m
clear all;
close all;

%% Initialize the toolbox
simulation = MotorProto('UQM145_SMPM');

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('SMPM Rotor','SynchronousRotor');
stator = model.newAssembly('SMPM Stator','Stator');

%% Define General Machine Parameters
nPoles            = 18;
nTeethPerPhase    = 2;
nTeeth            = 3 * nPoles * nTeethPerPhase;
nTurnsPerSlot     = 3;
len               = 0.1428877;
statorOuterRadius = 0.1250908;
statorInnerRadius = 0.0986028;
rotorOuterRadius  = 0.0971550;
rotorInnerRadius  = 0.0853906;
w_r               = 135;

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = len;
stator.Poles               = nPoles;
stator.Teeth               = nTeeth;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = Arnon7;
stator.SourceType          = 'CurrentSource';   %'CurrentSource', 'VoltageSource'
stator.CouplingType        = CouplingTypes.Dynamic;
stator.ConnectionType      = 'Wye';             %'Wye', 'Delta'
stator.WindingType         = 'Distributed';     %'Distributed', 'Concentrated'

%% Slot
stator.Slot.Turns = nTurnsPerSlot;

%% Stranded-Style Conductors
stator.Slot.ConductorType                 = 'Circular';
stator.Slot.Conductor.ConductorDiameter   = 1.0e-3 * 0.9 * 2 / (2^(0));
stator.Slot.Conductor.InsulationThickness = 0.10e-3 * 1.1 * 2 / (2^(0));

%% Bus-Bar Style Conductors
% stator.Slot.ConductorType           = 'Homogenized';
% stator.Slot.Conductor.PackingFactor = 0.5;

%% Define slot geometry
Bsat = 2;
Bmag = 1.23;

slotWidth   = 1-statorInnerRadius/rotorOuterRadius*Bmag/Bsat;
slotLength  = 1-statorInnerRadius*(1-slotWidth)*2*pi/nPoles*(3*nTeethPerPhase-1)/(3*nTeethPerPhase)/(statorOuterRadius-statorInnerRadius)/2;
notchWidth  = 0.1795;
notchLength = 0.01263;
[slotOutline, slotNotch] = slotTemplate(nTeeth, statorInnerRadius, statorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 1, 'auto', 'InnerSlotShape','rounded','OuterSlotShape','rounded');

stator.Slot.Shape        = slotOutline;
stator.ConductorMaterial = Copper;

stator.addRegion('slot', slotNotch, Air, DynamicsTypes.Static); 

%% Set Rotor Parameters
rotor.Poles               = nPoles;
rotor.Length              = len;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = Arnon7;
rotor.OperatingMode       = OperatingModes.Synchronous;
rotor.InitialAngle        = pi/nPoles*0;
rotor.BackironType        = BackironTypes.Laminated;

%% Create Rotor Permanent Magnet
pmRing     = 0.775e-3;
pmWidth    = 4e-3+pmRing;
pmEmbrace  = (1-0.11);
pmLength   = 2 * (rotorOuterRadius - pmWidth) * tan(2 * pi / nPoles / 2) * pmEmbrace;
pmPosition = [statorOuterRadius / 2 + rotorOuterRadius - pmWidth / 2, 0];

pmBody = Geometry2D.draw('Rect', 'Width', pmLength, 'Length', pmWidth + statorOuterRadius, 'Base', 'Center', 'Position', pmPosition, 'PlotStyle', {'m'});
pmTrim = Geometry2D.draw('Sector', 'Radius', [rotorInnerRadius, rotorOuterRadius-pmRing], 'Angle', 2 * pi / nPoles, 'Rotation', - pi / nPoles);
permanentMagnet = pmBody * pmTrim;

rotor.addRegion('pm', permanentMagnet,  NdFe35, DynamicsTypes.Dynamic);

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

rotor.addRegion('trimUHP', trimUHP, Air, DynamicsTypes.Static);
rotor.addRegion('trimLHP', trimLHP, Air, DynamicsTypes.Static);

retainingRing = Geometry2D.draw('Sector', 'Radius', [rotorOuterRadius-pmRing, rotorOuterRadius], 'Angle', 2 * pi / nPoles, 'Rotation', - pi / nPoles,'PlotStyle',{'b'});
rotor.addRegion('ring', retainingRing, Arnon7, DynamicsTypes.Static);

%% Set mesh parameters
mesh                       = simulation.Mesh;
mesh(1).MaximumElementSize = pmWidth / 4;
mesh(2).MaximumElementSize = (2*pi*statorInnerRadius)*(0.5/nTeeth)*0.28;

%% Set Excitation
% stator.Sources.ElectricalFrequency = w_r * nPoles / 2;

%% Voltage Source
%% %% 200 N-m, Maximum Current/Field Weakening
% stator.SourceType = 'VoltageSource';
% stator.Circuits.ElectricalFrequency = w_r * nPoles / 2;
% stator.Circuits.HarmonicNumbers     = 1;
% stator.Circuits.HarmonicAmplitudes  = 2*257.6445;
% stator.Circuits.HarmonicPhases      = -1.2727;

%% Current Source, 200 N-m, Maximum Current/Field Weakening
stator.SourceType = 'CurrentSource';
stator.Circuits.ElectricalFrequency = w_r * nPoles / 2;
stator.Circuits.HarmonicNumbers    = 1;
stator.Circuits.HarmonicAmplitudes = 250*1;
stator.Circuits.HarmonicPhases     = -2*pi/3 + pi*(1/4+1/8+1/16-1/32-1/64+1/128-1/256-1/512) + pi/2*0;

% stator.SourceType = 'DQOCurrentRegulator';
% stator.Circuits.ElectricalFrequency = w_r * nPoles / 2;
% stator.Circuits.InitialAngle        = -pi / 3;
% stator.Circuits.Idq = [101.5, -289];
% stator.Circuits.Kp  = [5, 5]*1e-4 * 75;
% stator.Circuits.Ki  = [5, 5]*1e-2 * 75;

model.build;
mesh.build;

nTimePoints = 47;
%simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'Verbose', true);
simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StorageLevel', 3, 'Verbose', true, 'MaxGMRESIterations', 25, 'ShootingTolerance', 1e-6, 'NewtonTolerance', 1e-6, 'GMRESTolerance', 1e-3, 'SymmetricJacobian', true);
%simulation.configureAlgorithm('TPFEM', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StorageLevel', 3, 'Verbose', true, 'MaxGMRESIterations', 5, 'NewtonTolerance', 1e-6, 'GMRESTolerance', 1e-3, 'SymmetricJacobian', true);

i = 1;
%nTimePoints =[  6,18,30,42,54
nTimePoints = [1,7,19,31,43];
simulation.configureAlgorithm('HarmonicBalance', 'TimePoints', nTimePoints(i));
solution = simulation.run;

%% Plotting
% solution.plot('A','Time',1);
% solution.plot('B','Time',1);
% solution.plot('H','Time',1);
% solution.plot('M','Time',1);
% solution.plot('A','Harmonic',[0, 1]);
% solution.plot('B','Harmonic',[0, 1]);
% solution.plot('H','Harmonic',[0, 1]);
% solution.plot('M','Harmonic',[0, 1]);
% solution.plot('LossDensity', 'UseSinglePlot', true);
solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% solution.plot('J','Time',1);
% solution.plot('J','Harmonic',1);
% solution.plot('E','Time',1);
% solution.plot('E','Harmonic',1);
% 
solution.plot('Flux Linkage','Time');
solution.plot('Flux Linkage','Harmonic');
solution.plot('Torque','Time');
solution.plot('Torque','Harmonic');
solution.plot('Voltage','Time');
solution.plot('Voltage','Harmonic');
solution.plot('Current','Time');
solution.plot('Current','Harmonic');