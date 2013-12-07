clear all;
close all;

display('IM_ROTOR_TEST_SCRIPT');
tic
%% Initialize the toolbox
RootObject = MotorProto('Induction Machien Rotor Tutorial');

%% Test adding a StatorLamination object to the Model
Model = RootObject.Model;
Tooth = Model.newComponent('myTooth','ToothComponent');

assert(Tooth == Model.Components(1));
assert(strcmp(class(Model.Components(1)),'ToothComponent'));
assert(strcmp(Model.Components(1).Name,'myTooth'));

%% Test removing a Pole
Model.removeComponent('myTooth')
assert(isempty(Model.Components));
assert(strcmp(Tooth.Name,'myTooth'));

%% Try adding two components with the same name
T1 = Model.newComponent('myTooth','ToothComponent');
try
	Model.newComponent('myTooth','ToothComponent');
end
assert(numel(Model.Components)==1);
assert(strcmp(Model.Components(1).Name,'myTooth'));

%% Try adding a component with a different name
T2 = Model.newComponent('myOtherTooth','ToothComponent');

%% Remove one component
assert(isvalid(T1));
assert(numel(Model.Components) == 2);
assert(strcmp(Model.Components(2).Name,'myOtherTooth'));

%% Try removing a component that's not there
Model.removeComponent('myTooth_1');
assert(isvalid(T1));
assert(isvalid(T2));
assert(numel(Model.Components) == 2);
assert(strcmp(Model.Components(1).Name,'myTooth'));
assert(strcmp(Model.Components(2).Name,'myOtherTooth'));

%% Try deleting components
Model.deleteComponent('myOtherTooth');
Model.removeComponent('myTooth');
assert( isvalid(T1));
assert(~isvalid(T2));
assert(isempty(Model.Components));

%% Try creating a new component with a non-string name
try
	T3 = Model.newComponent(1,'ToothComponent');
catch ME
    T3 = [];
end
assert(isempty(T3));
assert(numel(Model.Components) == 0);

T2 = Model.newComponent('1','ToothComponent');
assert(numel(Model.Components) == 1);
assert(strcmp(Model.Components(1).Name,'1'));

%% Try adding a component which has a null name
assert(numel(Model.Components) == 1);
Model.deleteComponent('1');
assert(~isvalid(T2));
assert(numel(Model.Components) == 0);

Tooth = Model.newComponent('myTooth','ToothComponent');

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general tooth parameters
Parameters.new('nTeeth',      96 );
Parameters.new('outerRadius', 0.085 );
Parameters.new('innerRadius', 0.05 );

Tooth.Teeth        = nTeeth;
Tooth.InnerRadius  = innerRadius;
Tooth.OuterRadius  = outerRadius;

assert(Tooth.SolutionSpatialFrequency  == 96);
assert(Tooth.GeometryFrequency         == 96);
assert(Tooth.SolutionHalfWaveSymmetry  == false);
assert(Tooth.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency     == 96);
assert(Model.HasHalfWaveSymmetry     == false);

Tooth.build;
assert(Tooth.SolutionSpatialFrequency  == 96);
assert(Tooth.GeometryFrequency         == 96);
assert(Tooth.SolutionHalfWaveSymmetry  == false);
assert(Tooth.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency     == 96);
assert(Model.HasHalfWaveSymmetry     == false);

%% Try assigning some pole properties
Tooth.DefaultMaterial = IronExampleMaterial;
try
    Tooth.DefaultMaterial = 1;
end
assert(isa(Tooth.DefaultMaterial,'IronExampleMaterial'));

Tooth.build;
assert(numel(Tooth.Regions) == 1);
assert(strcmp(Tooth.Regions.Name,'myTooth_DefaultRegion'));
assert(isa(Tooth.Regions.Geometry,'Sector'));
assert(all(Tooth.Regions.Geometry.Radius.Value == [0.05 0.085]));
assert(abs(Tooth.Regions.Geometry.Angle.Value - 2*pi/96) < eps*2*pi);
assert(abs(Tooth.Regions.Geometry.vRotation+pi/96) < eps*2*pi);
assert(isa(Tooth.Regions.Material,'IronExampleMaterial'));


%% define a symmetric slot shape
Parameters.new('slotInsetPercent',   0.02);
Parameters.new('slotWidthPercent',   0.5 );
Parameters.new('slotLengthPercent',  0.2 );
Parameters.new('toothGapPercent',    0.15);

toothWidthAngle     = 2 * pi / nTeeth;
slotOuterRadius     =  outerRadius * (1 - slotInsetPercent)...
                      +innerRadius * slotInsetPercent;

slotOuterPosition	= [slotOuterRadius, 0];
slotLength          = (outerRadius - innerRadius) * slotLengthPercent;
slotInnerRadius     = slotOuterRadius - slotLength;
slotOuterWidth     	= slotInnerRadius * sin(toothWidthAngle * slotWidthPercent);

toothWidthAngle     = 0;

slotPoints       = [ slotOuterPosition(1),...
                        slotOuterPosition(2) - slotOuterWidth / 2 ;
                        slotOuterPosition(1),...
                        slotOuterPosition(2) + slotOuterWidth / 2 ];
slotPoints       = [ slotPoints;
                     slotPoints(2,1) - slotLength * cos(toothWidthAngle / 2),...
                     slotPoints(2,2) - slotLength * sin(toothWidthAngle / 2);
                     slotPoints(1,1) - slotLength * cos(toothWidthAngle / 2),...
                     slotPoints(1,2) + slotLength * sin(toothWidthAngle / 2)];
                   
conductorOutline = Geometry2D.draw(  'Polygon2D',...
                                        'Points',   slotPoints,...
                                        'PlotStyle',{'y'});

distanceSlotBack	= (slotPoints(3,2) - slotPoints(4,2)) / 2;
adjustmentFactor    = cos(toothWidthAngle / 2)...
                     +tan(toothWidthAngle / 2) * (1 + sin(toothWidthAngle/2));
roundingInnerRadius	= distanceSlotBack / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotBack^2);
roundingPosition	= [slotPoints(3,1) + roundingInnerRadius, 0];
roundingAngle       = pi + toothWidthAngle;
roundingRotation	= (pi + toothWidthAngle) / 2;

roundingSector      = Geometry2D.draw('Sector',   ...
                                        'Radius',   [roundingInnerRadius,...
                                                        roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);
                                
conductorOutline	=   conductorOutline - roundingSector;

distanceSlotFront	= slotOuterWidth / 2;
adjustmentFactor    =  cos(toothWidthAngle/2)...
                      -tan(toothWidthAngle/2) * (1 - sin(toothWidthAngle/2));
roundingInnerRadius	= distanceSlotFront / adjustmentFactor;
roundingOuterRadius	= 1.1 * sqrt(roundingInnerRadius^2 + distanceSlotFront^2);
roundingPosition	= [ slotPoints(1,1) - roundingInnerRadius, 0 ];
roundingAngle       = pi - toothWidthAngle;
roundingRotation	= -( pi + toothWidthAngle ) / 2;

roundingSector      = Geometry2D.draw(	'Sector',   ...
                                        'Radius',   [roundingInnerRadius,...
                                                        roundingOuterRadius],...
                                        'Rotation', roundingRotation,...
                                        'Angle',    roundingAngle,...
                                        'Position', roundingPosition);

conductorOutline	= conductorOutline-roundingSector;

Tooth.addRegion( 'slot',...
                 conductorOutline,...
                 CopperExampleMaterial,...
                 'Dynamic',...
                 'Voltage');
Tooth.build;

assert(numel(Tooth.Regions) == 2);
assert(strcmp(Tooth.Regions(1).Name,'slot'));
assert(strcmp(Tooth.Regions(2).Name,'myTooth_DefaultRegion'));
assert(numel(Tooth.Components) == 2);
assert(strcmp(Tooth.Components(1).Name,'slot'));
assert(strcmp(Tooth.Components(2).Name,'myTooth_DefaultRegion'));
areaDifference = abs(    pi*(0.085^2-0.05^2)/96 ...
                        -Tooth.Regions(1).Geometry.area ...
                        -Tooth.Regions(2).Geometry.area);
assert(areaDifference < sqrt(eps)*2*pi/960);
assert(area(Tooth.Regions(1).Geometry*Tooth.Regions(2).Geometry) < sqrt(eps)*2*pi/80);
assert(Tooth.SolutionSpatialFrequency  == 96);
assert(Tooth.GeometryFrequency         == 96);
assert(Tooth.SolutionHalfWaveSymmetry  == false);
assert(Tooth.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency  == 96);
assert(Model.HasHalfWaveSymmetry       == false);

Tooth.removeComponent('slot');
assert(numel(Tooth.Regions) == 1);
assert(numel(Tooth.Components) == 1);
Tooth.build;

crookedSlot = rotate(copy(conductorOutline),pi/100,slotOuterPosition-slotLength/2);
Tooth.addRegion( 'cslot',...
                crookedSlot,...
                CopperExampleMaterial,...
                'Dynamic',...
                'Voltage');
Tooth.build;
            
assert(Tooth.SolutionSpatialFrequency  == 96);
assert(Tooth.GeometryFrequency         == 96);
assert(Tooth.SolutionHalfWaveSymmetry  == false);
assert(Tooth.GeometryMirrorSymmetry    == false);
assert(Model.SolutionSpatialFrequency     == 96);
assert(Model.HasHalfWaveSymmetry     == false);

Tooth.addRegion( 'slot',...
                 conductorOutline,...
                 CopperExampleMaterial,...
                 'Dynamic',...
                 'Voltage');
assert(numel(Tooth.Regions) == 3);
assert(numel(Tooth.Components) == 3);

Tooth.removeComponent( 'cslot');
assert(numel(Tooth.Regions) == 2);
assert(numel(Tooth.Components) == 2);

Parameters.new('fbMinBridgeWidth', 0.0004);

gapOuterRadius = outerRadius;
gapInnerRadius = outerRadius - roundingInnerRadius;
gapWidthAngle  = 2 * pi / nTeeth * toothGapPercent;

gap            = Geometry2D.draw(	'Sector',   ...
                                 	'Radius',	[gapInnerRadius,...
                                                        gapOuterRadius],...
                                  	'Angle',  	gapWidthAngle,...
                                   	'Rotation',	-gapWidthAngle/2,...
                                  	'PlotStyle',{'w'});

gap            = gap - conductorOutline;
null           = gap * conductorOutline;
assert(isempty(null.Curves));
assert(null.area == 0);

Tooth.addRegion( 'tgap',...
                gap,...
                Air,...
                'Static',...
                'Current');
assert(numel(Tooth.Regions) == 3);
Tooth.build;

rotGap  = rotate(copy(gap),2*pi*95/96,[0 0]);
rotSlot = rotate(copy(conductorOutline),2*pi*95/96,[0 0]);
s       = rotGap + rotSlot;

assert(Tooth.SolutionSpatialFrequency  == 96);
assert(Tooth.GeometryFrequency         == 96);
assert(Tooth.SolutionHalfWaveSymmetry  == false);
assert(Tooth.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency  == 96);
assert(Model.HasHalfWaveSymmetry       == false);

for i = 1:8
    Tooth = Tooth.build(i/96);
    assert(numel(Tooth.Regions) == 2*i+1);
    assert(numel(Tooth.Components) == 2*i+1);
end

%% try building full model
build(Tooth,1);

%% Test Run
RootObject.configureAlgorithm('Static');
RootObject.Model.Mesh.MaximumElementSize = (outerRadius - innerRadius) / 20;
RootObject.Model.Mesh.MinimumElementSize = 0;
s = RootObject.run;

%% test various matrix properties
K = s.Matrices.K(0);
isSymmetric        = max(max(abs(K-K.'))) / max(max(abs(K))) < eps;
[~,p]              = chol(K);
isPositiveDefinite = (p == 0);
assert(isPositiveDefinite);
assert(isSymmetric);

%% Todo
%   In PoleComponent, the test for GeometryMirrorSymmetry does not take into
%   acount which regions are slave regions and which regions are not. This
%   should be changed immediately.
%
%   Begin development of new stator model, factoring some of the dependent
%   properties out of PoleComponent
%
%   We could seperate the construction of the geometry from the construction of
%   the windings and input functions, and specify an inferiorto relationship
%   between a PoleComponent and ToothComponent
%
%   We could specify that slave components are in series or in parallel
%
%   Think about extending PoleComponent to add specific addRegion methods
%   for different materials. E.g. addPermanentMagnet, addFluxBarrier
%
%   Write automated tests for overlapping components <- Component
%
%   Write tests for geometry which overlap on edges
%       -Write more geometry tests
toc 