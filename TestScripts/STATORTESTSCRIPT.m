clear all;close all;
tic
%% Initialize the toolbox
RootObject = MotorProto('Stator Model Tutorial');

%% Test adding a StatorLamination object to the Model
Stator = RootObject.newComponent('myStator','StatorLamination');

assert(Stator == RootObject.Model(1));
assert(strcmp(class(RootObject.Model(1)),'StatorLamination'));
assert(strcmp(RootObject.Model(1).Name,'myStator'));

%% Test removing a StatorLamination methods
assert(isempty(Stator.remove('myStator')));
assert(Stator.isNamed('myStator'));

%% Test removing a StatorLamination from the model
RootObject.removeComponent('myStator');
assert(isempty(RootObject.Model));

%% Try adding two components with the same name
S1 = RootObject.newComponent('myStator','StatorLamination');
try
    S2 = RootObject.newComponent('myStator','StatorLamination'); 
end
assert(numel(RootObject.Model)==1);

try
    RootObject.addComponent(S1); 
end
assert(numel(RootObject.Model)==1);

%% Try adding a component with a different name
S2 = RootObject.newComponent('myOtherStator','StatorLamination');

%% Remove one component
RootObject.removeComponent('myStator');
assert(isvalid(S1));
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myOtherStator'));

%% Try removing a component that's not there
RootObject.removeComponent('myStator');
assert(isvalid(S1));
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myOtherStator'));

%% Try deleting components
RootObject.deleteComponent('myOtherStator');
assert(~isvalid(S2));
assert(isempty(RootObject.Model));

%% Try adding components
RootObject.addComponent(S1);
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myStator'));

%% Try creating a new component with a non-string name
try
    S2 = RootObject.newComponent(1,'StatorLamination');
end
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myStator'));

S2 = RootObject.newComponent('1','StatorLamination');
assert(numel(RootObject.Model) == 2);
assert(RootObject.Model(2).isNamed('1'));

%% Try adding a component which has a null name
S3 = StatorLamination;
try
    RootObject.addComponent(S3);
end
assert(numel(RootObject.Model) == 2);
RootObject.deleteComponent('1');

Stator                 = RootObject.Model(1);

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general stator parameters
Parameters.new('nPoles',            8   );
Parameters.new('nTeethPerPhase',	3   );
Parameters.new('nTurns',            2   );
Parameters.new('outerRadius',       0.12);
Parameters.new('innerRadius',       0.08513);

%% define slot parameters
Parameters.new('slotInsetPercent',   0.02);
Parameters.new('slotWidthPercent',   0.5 );
Parameters.new('slotLengthPercent',  0.6 );
Parameters.new('toothGapPercent',    0.15 );

Stator.Poles        = nPoles;
Stator.Teeth        = 3 * nPoles * nTeethPerPhase;
Stator.Turns        = nTurns;
Stator.InnerRadius  = innerRadius;
Stator.OuterRadius  = outerRadius;

assert(Stator.SolutionSpatialFrequency  == 4);
assert(Stator.SolutionSpatialSymmetry   == 2);
assert(Stator.GeometricFrequency        == 72);

%% Try assigning some stator properties
Stator.DefaultMaterial = IronExampleMaterial;
Stator.SlotMaterial    = CopperExampleMaterial; 
try
    Stator.DefaultMaterial = 1;
end
assert(isa(Stator.DefaultMaterial,'MaterialProperty'));

try
    Stator.SlotMaterial = 1;
end
assert(isa(Stator.SlotMaterial,'MaterialProperty'));

%% define rough slot geometry
toothWidthAngle     = 2 * pi / 3 / nPoles / nTeethPerPhase;
slotInnerRadius     =  innerRadius * (1 - slotInsetPercent)...
                      +outerRadius * slotInsetPercent;
slotInnerPosition	= [slotInnerRadius, 0];
slotInnerWidth     	= slotInnerRadius * sin(toothWidthAngle * slotWidthPercent);

slotLength          = (outerRadius - innerRadius) * slotLengthPercent;

toothWidthAngle     = 0;

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
roundingPosition	= [ slotPoints(1,1) + roundingInnerRadius, 0 ];
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
try
    Stator.SlotOutline = 1;
end
assert(isa(Stator.SlotOutline,'Geometry2D'));

%% plot
% Stator.build;
% figure;
% Stator.plot;
% pause(1);

% %% setup simulation
% %% Current Drive, Static Simulation
% Stator.Input                = 'Current';
% Stator.Dynamics             = 'Static';
% Stator.FundamentalFrequency = 2*pi*60;                    
% Stator.HarmonicAmplitudes   = 200;
% Stator.HarmonicAngles      	= pi-2*pi/3;
% Stator.HarmonicNumbers      = 1;
% Stator.ExternalResistance   = 1e-3;
% Stator.ExternalInductance   = 1e-6;
% 
% %% assign mesh parameters
% Stator.MinimumElementQuality     = 1;
% Stator.MaximumBoundaryEdgeLength = 2*innerRadius*sin(pi/288);
% Stator.InternalEdgeUniformity    = 1.1;
% 
% Stator.generateMesh;
% figure;
% Stator.Mesh.plot;
% toc
% %pause;
% 
% % phase = Stator.build('phase');
% % clf;
% % n = numel(phase);
% % for i=1:n;
% % 	phase{i}.plot;
% % end
% % axis equal;
% % pause;
% 
% 
% % pole = Stator.build('pole');
% % clf;
% % n = numel(pole);
% % for i=1:n;
% %     pole{i}.plot;
% % end
% % axis equal;
% % P.edit('np',4);
toc 