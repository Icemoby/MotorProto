clear all;close all;clc;
%% Initialize the toolbox
RootObject = MotorProto('Stator Model Tutorial');

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general stator parameters
Parameters.new('nPoles',            6       );
Parameters.new('nTeethPerPhase',	3       );
Parameters.new('nTurns',            5       );
Parameters.new('outerRadius',       0.12    );
Parameters.new('innerRadius',       0.08513 );

%% define slot parameters
Parameters.new('slotInsetPercent',   0.02);
Parameters.new('slotWidthPercent',   0.4 );
Parameters.new('slotLengthPercent',  0.6 );
Parameters.new('toothGapPercent',    0.1 );

%% assign stator parameters;
Stator              = RootObject.Stator;
Stator.Poles        = nPoles;
Stator.Teeth        = 3 * nPoles * nTeethPerPhase;
Stator.Turns        = nTurns;
Stator.InnerRadius  = innerRadius;
Stator.OuterRadius  = outerRadius;

%% assign material properties
Stator.LaminationMaterial = LaminatedIronMaterial;
Stator.ConductorMaterial  = HomogenizedSlotCopper;

%% define rough slot geometry
toothWidthAngle     = 2 * pi / 3 / nPoles / nTeethPerPhase;
slotInnerRadius     =  innerRadius * (1 - slotInsetPercent)...
                      +outerRadius * slotInsetPercent;
slotInnerPosition	= [slotInnerRadius, 0];
slotInnerWidth     	= slotInnerRadius * sin(toothWidthAngle * slotWidthPercent);

slotLength          = (outerRadius - innerRadius) * slotLengthPercent;

%% polygon slot
slotPoints       = [ slotInnerPosition(1),...
                        slotInnerPosition(2) + slotInnerWidth / 2 ;
                        slotInnerPosition(1),...
                        slotInnerPosition(2) - slotInnerWidth / 2 ];
slotPoints       = [ slotPoints;
                     slotPoints(2,1) + slotLength * cos(toothWidthAngle / 2),...
                     slotPoints(2,2) - slotLength * sin(toothWidthAngle / 2);
                     slotPoints(1,1) + slotLength * cos(toothWidthAngle / 2),...
                     slotPoints(1,2) + slotLength * sin(toothWidthAngle / 2)];
                   
conductorOutline = Geometry2D.draw(  'Polygon2D',...
                                        'Points',   slotPoints,...
                                        'PlotStyle',{'y'});

%% round outer edges
distanceSlotBack	= (slotPoints(4,2) - slotPoints(3,2)) / 2;
adjustmentFactor    = cos(toothWidthAngle / 2)...
                     +tan(toothWidthAngle / 2) * (1 + sin(toothWidthAngle/2));
roundingInnerRadius	= distanceSlotBack / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotBack^2);
roundingPosition	= [slotPoints(3,1) - roundingInnerRadius, 0];
roundingAngle       = pi + toothWidthAngle;
roundingRotation	= -(pi + toothWidthAngle) / 2;

roundingSector      = Geometry2D.draw('Sector',   ...
                                        'Radius',   [roundingInnerRadius,...
                                                        roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);
                                
conductorOutline	=   conductorOutline - roundingSector;

%% round inner edges
distanceSlotFront	= slotInnerWidth / 2;
adjustmentFactor    =  cos(toothWidthAngle/2)...
                      -tan(toothWidthAngle/2) * (1 - sin(toothWidthAngle/2));
roundingInnerRadius	= distanceSlotFront / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotFront^2);
roundingPosition	= [slotPoints(1,1) + roundingInnerRadius, 0 ];
roundingAngle       = pi - toothWidthAngle;
roundingRotation	= ( pi + toothWidthAngle ) / 2;

roundingSector      = Geometry2D.draw(	'Sector',   ...
                                        'Radius',   [roundingInnerRadius,...
                                                        roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);

conductorOutline	= conductorOutline-roundingSector;

%% define gap between teeth
gapInnerRadius      = innerRadius;
gapOuterRadius      = innerRadius + roundingInnerRadius;
gapWidthAngle       = 2 * pi / 3 / nPoles / nTeethPerPhase * toothGapPercent;

gap                 = Geometry2D.draw(	'Sector',   ...
                                        'Radius',	[gapInnerRadius,...
                                                        gapOuterRadius],...
                                        'Angle',  	gapWidthAngle,...
                                        'Rotation',	-gapWidthAngle/2);

conductorOutline	= conductorOutline + gap;

%% save conductorOutline to stator
Stator.SlotOutline = conductorOutline;

figure;
plot(Stator);axis equal;
pause;

%% generate the stator mesh
Stator = generateMesh(Stator);

clf;
plot(Stator.Mesh);axis equal;
pause;

%% setup a solver

RootObject.Solver = Solvers.HarmonicBalanceDemo
RootObject.run;

%% plot a flux density snapshot
RootObject.Solution.plot('B','Time',0);

%% plot a few current density harmonic magnitudes to show high frequency effects
RootObject.Solution.plot('J','Harmonics',[1 3;5 7]);

%% check average torque and losses
RootObject.Solution.AverageTorque
RootObject.Solution.TotalLosses

%Stator.generateMesh;
% %pause;
% 
% clf;
% Stator.plot;
% %pause;
% 
% phase = Stator.build('phase');
% clf;
% n = numel(phase);
% for i=1:n;
% 	phase{i}.plot;
% end
% axis equal;
%pause;
% 
% 
% pole = Stator.build('pole');
% clf;
% n = numel(pole);
% for i=1:n;
%     pole{i}.plot;
% end
% axis equal;
% % P.edit('np',4);
% toc 