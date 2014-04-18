%% UQM145_SMPM_Tutorial.m
clear all;
close all;

tic

%% Initialize the toolbox
simulation = MotorProto('UQM145_SMPM');
%warning on 'MotorProto:Verbose';

%% Add components to the model
model  = simulation.Model;
rotor  = model.newAssembly('SMPM Rotor','SynchronousRotor');
stator = model.newAssembly('SMPM Stator','Stator');

%% Define General Machine Parameters
nPoles            = 18;
nTeethPerPhase    = 2;
nTeeth            = 3 * nPoles * nTeethPerPhase;
len               = 0.143;
statorOuterRadius = 0.125;
statorInnerRadius = 0.0986;
rotorOuterRadius  = 0.0972;
rotorInnerRadius  = 0.0854;
w_r               = 115;

%% Define Stator Geometry and Material Properties
stator.ElectricalFrequency = w_r * nPoles / 2;
stator.Length              = len;
stator.Poles               = nPoles;
stator.Teeth               = nTeeth;
stator.InnerRadius         = statorInnerRadius;
stator.OuterRadius         = statorOuterRadius;
stator.DefaultMaterial     = Arnon7;
stator.SourceType          = 'VoltageSource';
stator.ConnectionType      = 'Wye';
stator.ConductorDynamics   = 'Dynamic';
stator.WindingType         = 'Distributed';
stator.Slot.Turns          = 3;

%%Stranded-Style Conductors
stator.Slot.ConductorType                 = 'Circular';
stator.Slot.Conductor.ConductorDiameter   = 2.1e-3 / 1;
stator.Slot.Conductor.InsulationThickness = 0.21e-3 / 1;

%%Bus-Bar Style Conductors
% stator.Slot.ConductorType           = 'Homogenized';
% stator.Slot.Conductor.PackingFactor = 0.5;

%% Define slot geometry
slotWidth   = 0.5;
slotLength  = 0.69;
notchWidth  = (1-slotWidth)  * 0.4;
notchLength = (1-slotLength) * 0.15;

[slotOutline, slotNotch] = slotTemplate(nTeeth, statorInnerRadius, statorOuterRadius, notchWidth, notchLength, slotWidth, slotLength, 1, 'auto','InnerSlotShape','rounded','OuterSlotShape','rounded');

stator.Slot.Shape        = slotOutline;
stator.ConductorMaterial = CopperExampleMaterial;

stator.addRegion('slot', slotNotch, Air, 'Static', -1); 

%% Set Rotor Parameters
rotor.Poles               = nPoles;
rotor.Length              = len;
rotor.ElectricalFrequency = w_r * nPoles / 2;
rotor.InnerRadius         = rotorInnerRadius;
rotor.OuterRadius         = rotorOuterRadius;
rotor.DefaultMaterial     = Arnon7;
rotor.OperatingMode       = 'synchronous';
rotor.InitialAngle        = 0;

%% Create Rotor Permanent Magnet
pmRing     = 0.775e-3;
pmWidth    = 4e-3+pmRing;
pmEmbrace  = (1-0.11);
pmLength   = 2 * (rotorOuterRadius - pmWidth) * tan(2 * pi / nPoles / 2) * pmEmbrace;
pmPosition = [statorOuterRadius / 2 + rotorOuterRadius - pmWidth / 2, 0];

pmBody = Geometry2D.draw('Rect', 'Width', pmLength, 'Length', pmWidth + statorOuterRadius, 'Base', 'Center', 'Position', pmPosition, 'PlotStyle', {'m'});
pmTrim = Geometry2D.draw('Sector', 'Radius', [rotorInnerRadius, rotorOuterRadius-pmRing], 'Angle', 2 * pi / nPoles, 'Rotation', - pi / nPoles);
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
mesh(1).MaximumElementSize = pmWidth / 4;
mesh(2).MaximumElementSize = (2*pi*statorInnerRadius)*(0.5/nTeeth)*0.28;

%% Set Excitation
stator.Sources.ElectricalFrequency = w_r * nPoles / 2;

%% Voltage Source
%% w_r = 115Hz, Max Current, 145kW/Min Voltage
stator.SourceType = 'VoltageSource';
stator.Sources.HarmonicNumbers    = 1;
stator.Sources.HarmonicAmplitudes = 466;
stator.Sources.HarmonicPhases     = (180+105)*(pi/180);

%% Max Current/Max Torque, 145kW/w_r = 3500/60
% stator.SourceType = 'VoltageSource';
% stator.Sources.HarmonicNumbers    = 1;
% stator.Sources.HarmonicAmplitudes = 410;
% stator.Sources.HarmonicPhases     = (180+110)*(pi/180);

%% Current Source
%% w_r = 115Hz, Max Current, 145kW/Min Voltage
% stator.SourceType = 'CurrentSource';
% i = 145 * exp(1i*(7/6)*pi) * 0.36;
% i = i + 145 * exp(1i*10/6*pi) * sqrt(1-0.36^2) * 1;
% stator.Sources.HarmonicNumbers    = 1;
% stator.Sources.HarmonicAmplitudes = abs(i)*1;
% stator.Sources.HarmonicPhases     = angle(i);

%% Max Current/Max Torque, 145kW/w_r = 3500/60
% stator.SourceType = 'CurrentSource';
% i = 145 * exp(1i*(7/6)*pi) * 0.77;
% i = i + 145 * exp(1i*10/6*pi) * sqrt(1-0.77^2);
% stator.Sources.HarmonicNumbers    = 1;
% stator.Sources.HarmonicAmplitudes = abs(i);
% stator.Sources.HarmonicPhases     = angle(i);

% simulation.configureAlgorithm('Static', 'TimePoints', nTimePoints, 'Verbose', true);
% solution = simulation.run;%         

% nTimePoints = 2^1;
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 4, 'StoreDecompositions', true, 'Verbose', true, 'TransientSolver', false, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-6);
% solution = simulation.run;
% nTimePoints = 2^9;
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', 4, 'StoreDecompositions', true, 'Verbose', true, 'TransientSolver', false, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-6);
% solution = simulation.run(solution.Algorithm.X(:,1));

% simulation.configureAlgorithm('TPFEM', 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', true, 'Verbose', true, 'MaxGMRESIterations', 5, 'NewtonTolerance', 1e-3);
% solution = simulation.run;

% N           = 7;
% nStages    	= [1 2 4];
% M           = length(nStages);
% losses      = zeros(N,M,3);
% a_0        	= cell(N,M,3);
% a_1      	= cell(N,M,3);
% x_0        	= cell(N,M,3);
% sim_time  	= zeros(N,M,3);
% n          	= zeros(N,1);
% 
% for i = 1:N
%     nTimePoints = 2^(i+3);
%     for j = 1:M
%         %% Shooting Newton
%         k = 1;
%         simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', nStages(j), 'StoreDecompositions', true, 'Verbose', true, 'TransientSolver', false, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-6);
%         solution = simulation.run;
%         
%         l = solution.getBulkVariableData('AverageLosses');
%         losses(i,j,k) = l{1}(1);
%         
%         a = solution.getContinuumVariableData('A','Harmonic',0);
%         a_0{i,j,k} = a{1};
%         
%         a = solution.getContinuumVariableData('A','Harmonic',1);
%         a_1{i,j,k} = a{2};
%         
%         x_0{i,j,k} = solution.Algorithm.X(:,1);
%         
%         sim_time(i,j,k) = solution.Algorithm.SimulationTime;
%         
%         n(i) = length(solution.Algorithm.Times)-1;
%         
%         %% TPFEM
%         k = 2;
%         simulation.configureAlgorithm('TPFEM', 'TimePoints', nTimePoints, 'RungeKuttaStages', nStages(j), 'StoreDecompositions', true, 'Verbose', true, 'MaxGMRESIterations', 5, 'NewtonTolerance', 1e-6);
%         solution = simulation.run;
%         
%         l = solution.getBulkVariableData('AverageLosses');
%         losses(i,j,k) = l{1}(1);
%         
%         a = solution.getContinuumVariableData('A','Harmonic',0);
%         a_0{i,j,k} = a{1};
%         
%         a = solution.getContinuumVariableData('A','Harmonic',1);
%         a_1{i,j,k} = a{2};
%         
%         x_0{i,j,k} = solution.Algorithm.X(:,1);
%         
%         sim_time(i,j,k) = solution.Algorithm.SimulationTime;
% 
%         %% Transient
%         k = 3;
%         simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', nStages(j), 'StoreDecompositions', false, 'Verbose', true, 'TransientSolver', true, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-4);
%         solution = simulation.run;
%         
%         l = solution.getBulkVariableData('AverageLosses');
%         losses(i,j,k) = l{1}(1);
%         
%         a = solution.getContinuumVariableData('A','Harmonic',0);
%         a_0{i,j,k} = a{1};
%         
%         a = solution.getContinuumVariableData('A','Harmonic',1);
%         a_1{i,j,k} = a{2};
%         
%         x_0{i,j,k} = solution.Algorithm.X(:,1);
%         
%         sim_time(i,j,k) = solution.Algorithm.SimulationTime;
%         
%         %% Save
%         save('C:\\results','losses','a_0','a_1','x_0','sim_time');
%         pause(1);
%     end
% end
% 
% %% Reference Solution
% nTimePoints = 2^((3+N)+4);
% simulation.configureAlgorithm('ShootingNewton', 'TimePoints', nTimePoints, 'RungeKuttaStages', nStages(M), 'StoreDecompositions', false, 'Verbose', true, 'TransientSolver', false, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-6);
% solution = simulation.run(x_0{end,end,1});
% 
% l = solution.getBulkVariableData('AverageLosses');
% ref_losses = l{1}(1);
% 
% a = solution.getContinuumVariableData('A','Harmonic',0);
% ref_a_0 = a{1};
% 
% a = solution.getContinuumVariableData('A','Harmonic',1);
% ref_a_1 = a{2};
% 
% ref_x_0 = solution.Algorithm.X(:,1);
% 
% ref_sim_time = solution.Algorithm.SimulationTime;
% 
% save('C:\\results','losses','a_0','a_1','x_0','sim_time','ref_losses','ref_a_0','ref_a_1','ref_x_0','ref_sim_time');

% simulation.configureAlgorithm('ShootingNewton', 'TransientSolver', true, 'TimePoints', nTimePoints, 'RungeKuttaStages', 2, 'StoreDecompositions', false, 'Verbose', true, 'MaxGMRESIterations', 5, 'ShootingTolerance', 1e-5);
% solution = simulation.run;

% simulation.configureAlgorithm('HarmonicBalanceDomainDecomposition', 'TimePoints', nTimePoints, 'Verbose', true, 'NewtonTolerance', 1e-2, 'GMRESTolerance', 1e-2);
% solution = simulation.run;

%l = solution.getBulkVariableData('AverageConductionLosses');
%t = solution.getBulkVariableData('Torque','Harmonic');

%% Plotting
% solution.plot('A','Time',1);
% solution.plot('B','Time',1);
% % solution.plot('H','Time',1);
% % solution.plot('M','Time',1);
% % solution.plot('A','Harmonic',[0, 1]);
% solution.plot('B','Harmonic',[0, 1]);
% solution.plot('LossDensity', 'UseSinglePlot', true, 'DataFunction', @(x)(log10(x)), 'DataFunctionString', 'log_{10}');
% % % % solution.plot('J','Harmonic',1);
% % % % solution.plot('J','Time',1);
% % % % 
% % % solution.plot('Flux Linkage','Time');
% % % solution.plot('Flux Linkage','Harmonic');
% % % solution.plot('Torque','Time');
% solution.plot('Torque','Harmonic');
% solution.plot('Voltage','Time');
% solution.plot('Voltage','Harmonic');
% solution.plot('Current','Time');
% solution.plot('Current','Harmonic');