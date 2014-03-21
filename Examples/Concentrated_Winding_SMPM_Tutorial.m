%%% Call this part once %%%
%% Concentrated_Winding_SMPM_Tutorial.m
% clear all;
% close all;

%% Initialize the toolbox
simulation = MotorProto('Concentrated Winding SMPM');

%% Add components to the model
model  = simulation.Model;

%%% Put the rest within an optimizaiton loop %%%
%% Define General Machine Parameters
f_r               = 70;
nPoles            = 50;
nTeeth            = 60;

stackLength       = 0.14; %abs

statorInnerRadius = 0.1405; %abs
statorBILength    = 0.01;   %abs
statorBICutRadius = 0.005;  %abs

toothYokeLength   = 0.02; %abs
toothYokeWidth    = 0.5;  %pct theta
toothFaceLength   = 0.003;%abs
toothFaceChamfer  = 0.5;  %pct toothFaceLength
toothGapWidth     = 0.01; %pct theta
turnsPerTooth     = 1;    %abs

slotFilletSize1   = 0.5;  %pct radius
slotFilletSize2   = 1;    %pct radius
slotPadding       = 0.001;%abs
slotPackingFactor = 0.5;  %pct conductor

airgapLength      = 0.001;%abs

pmLength          = 4e-3; %abs
pmEmbrace         = 0.94; %pct theta
pmFilletSize      = 0.5;  %pct radius

rotorBILength     = 0.013;%abs
rotorBICutLength  = 0.007;%abs
rotorBICutWidth   = 0.9;  %pct

statorIronMaterial = Steel1010;
rotorIronMaterial  = Steel1010;
pmMaterial         = NdFe35;

[model, stator, rotor] = make_CW_SMPM_Machine(model,f_r,nPoles,nTeeth,stackLength,turnsPerTooth,statorInnerRadius, statorBILength, statorBICutRadius, toothYokeLength, toothYokeWidth, toothFaceLength, toothFaceChamfer, toothGapWidth, slotFilletSize1, slotFilletSize2, slotPadding, slotPackingFactor,airgapLength, pmLength, pmEmbrace, pmFilletSize,rotorBILength, rotorBICutWidth, rotorBICutLength, statorIronMaterial, rotorIronMaterial, pmMaterial);

%% Set mesh parameters
mesh = simulation.Mesh;
mesh(1).MaximumElementSize = statorBILength / 3;
mesh(2).MaximumElementSize = rotorBILength / 3;

%% Voltage Source
% stator.SourceType = 'VoltageSource';
% stator.Sources.HarmonicNumbers    = 1;
% stator.Sources.HarmonicAmplitudes = 708;
% stator.Sources.HarmonicPhases     = -2 * pi * 101 / 360;

%% Current Source
stator.SourceType                 = 'CurrentSource';
stator.Sources.HarmonicNumbers    = 1;
stator.Sources.HarmonicAmplitudes = 225 / sqrt(3);
stator.Sources.HarmonicPhases     = 0;

%% Configure algorithm
% timePointsPerPeriod = 10;
% simulation.configureAlgorithm('Static', 'TimePoints', timePointsPerPeriod, 'Verbose', true);
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true,'ShootingTolerance',1e-4);
% solution = simulation.run;

%% Plotting
% solution.plot('A','Time',1);
% solution.plot('B','Time',1);
% solution.plot('H','Time',1);
% solution.plot('H','Harmonic',0);
% solution.plot('A','Harmonic',[0, model.TemporalSubharmonics]);
% solution.plot('B','Harmonic',[0, model.TemporalSubharmonics]);
% solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% solution.plot('J','Harmonic',model.TemporalSubharmonics);
% solution.plot('J','Time',1);
% solution.plot('E','Time',1);
% 
% solution.plot('FluxLinkage','Time');
% solution.plot('FluxLinkage','Harmonic');
% solution.plot('Torque','Time');
% solution.plot('Torque','Harmonic');
% solution.plot('Voltage','Time');
% solution.plot('Voltage','Harmonic');
% solution.plot('Current','Time');
% solution.plot('Current','Harmonic');

%% Data
% solution.getContinuumVariableData('H','Time',1)
% H      = getPMFieldIntensity(solution,pmMaterial);

% torque = solution.getBulkVariableData('Torque','Time');
% torque = torque{1};
%torque = solution.getBulkVariableData('Torque','Harmonic');

% flux_linkage = solution.getBulkVariableData('FluxLinkage','Time');
% flux_linkage = flux_linkage{1}{1};
%flux_linkage = solution.getBulkVariableData('Flux Linkage','Harmonic');

% mass       = solution.Model.Mass;
% statorMass = solution.Model.Assemblies(1).Mass;
% rotorMass  = solution.Model.Assemblies(2).Mass;

% model.build;
% model.plot;