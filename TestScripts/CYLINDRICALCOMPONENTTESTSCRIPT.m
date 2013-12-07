clear all;
close all;
clc;

tic
%% Initialize the toolbox
RootObject = MotorProto('Stator Model Tutorial');

%% Test adding a StatorLamination object to the Model
Pole = RootObject.newComponent('myPole','PoleComponent');

assert(Pole == RootObject.Model(1));
assert(strcmp(class(RootObject.Model(1)),'PoleComponent'));
assert(strcmp(RootObject.Model(1).Name,'myPole'));

%% Test removing a StatorLamination methods
assert(isempty(Pole.remove('myPole')));
assert(Pole.isNamed('myPole'));

%% Test removing a StatorLamination from the model
RootObject.removeComponent('myPole');
assert(isempty(RootObject.Model));

%% Try adding two components with the same name
P1 = RootObject.newComponent('myPole','PoleComponent');
try
    P2 = RootObject.newComponent('myPole','PoleComponent'); 
end
assert(numel(RootObject.Model)==1);

try
    RootObject.addComponent(P1); 
end
assert(numel(RootObject.Model)==1);

%% Try adding a component with a different name
P2 = RootObject.newComponent('myOtherPole','PoleComponent');

%% Remove one component
RootObject.removeComponent('myPole');
assert(isvalid(P1));
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myOtherPole'));

%% Try removing a component that's not there
RootObject.removeComponent('myPole');
assert(isvalid(P1));
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myOtherPole'));

%% Try deleting components
RootObject.deleteComponent('myOtherPole');
assert(~isvalid(P2));
assert(isempty(RootObject.Model));

%% Try adding components
RootObject.addComponent(P1);
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myPole'));

%% Try creating a new component with a non-string name
try
    P2 = RootObject.newComponent(1,'PoleComponent');
end
assert(numel(RootObject.Model) == 1);
assert(RootObject.Model.isNamed('myPole'));

P2 = RootObject.newComponent('1','PoleComponent');
assert(numel(RootObject.Model) == 2);
assert(RootObject.Model(2).isNamed('1'));

%% Try adding a component which has a null name
P3 = PoleComponent;
try
    RootObject.addComponent(P3);
end
assert(numel(RootObject.Model) == 2);
RootObject.deleteComponent('1');

Pole = RootObject.Model(1);

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general pole parameters
Parameters.new('nPoles',            8   );
Parameters.new('outerRadius',       0.085);
Parameters.new('innerRadius',       0.05);

Pole.Poles        = nPoles;
Pole.InnerRadius  = innerRadius;
Pole.OuterRadius  = outerRadius;

assert(Pole.SolutionSpatialFrequency  == 4);
assert(Pole.GeometryFrequency         == 8);
assert(Pole.SolutionHalfWaveSymmetry  == true);
assert(Pole.GeometryMirrorSymmetry    == true);
assert(Pole.ModelSpatialFrequency     == 4);
assert(Pole.ModelHalfWaveSymmetry     == true);

Pole.compile;
assert(Pole.SolutionSpatialFrequency  == 4);
assert(Pole.GeometryFrequency         == 8);
assert(Pole.SolutionHalfWaveSymmetry  == true);
assert(Pole.GeometryMirrorSymmetry    == true);
assert(Pole.ModelSpatialFrequency     == 4);
assert(Pole.ModelHalfWaveSymmetry     == true);

%% Try assigning some pole properties
Pole.DefaultMaterial = IronExampleMaterial;
try
    Pole.DefaultMaterial = 1;
end
assert(isa(Pole.DefaultMaterial,'IronExampleMaterial'));

Pole.compile;
assert(numel(Pole.Regions) == 1);
assert(strcmp(Pole.Regions.Name,'DefaultRegion'));
assert(isa(Pole.Regions.Geometry,'Sector'));
assert(all(Pole.Regions.Geometry.Radius.Value == [0.05 0.085]));
assert(abs(Pole.Regions.Geometry.Angle.Value - 2*pi/8) < eps*2*pi);
assert(Pole.Regions.Geometry.Rotation.Value == -pi/8);
assert(isa(Pole.Regions.Material,'IronExampleMaterial'));

%% Add a symmetric shape to the pole

%% define slot parameters
Parameters.new('nTeethPerPhase',	3   );
Parameters.new('nTurns',            2   );
Parameters.new('slotInsetPercent',   0.02);
Parameters.new('slotWidthPercent',   0.5 );
Parameters.new('slotLengthPercent',  0.6 );
Parameters.new('toothGapPercent',    0.15 );

%% define rough slot geometry
Parameters.new('pmInsetPercent',    0.00);
Parameters.new('pmLengthPercent',	0.83);
Parameters.new('pmWidthPercent',    0.40);

pmOuterRadius    = (outerRadius*cos(pi/nPoles))*(1-pmInsetPercent)+innerRadius*pmInsetPercent;
pmMaxLength      = 2*pmOuterRadius*tan(pi/nPoles);
pmLength         = pmMaxLength*pmLengthPercent;
pmMaxWidth       = pmOuterRadius - pmLength / ( 2 * tan( pi / nPoles ) );
pmWidth          = pmMaxWidth*pmWidthPercent;
pmPosition       = [pmOuterRadius-pmWidth/2 , 0];

permanentMagnet  = Geometry2D.draw(  'Rect',...
                                        'Width',    pmLength,...
                                        'Length',	pmWidth,...
                                        'Base',     'Center',...
                                        'Position', pmPosition,...
                                        'PlotStyle',{'m'});

Pole.addRegion( 'pm',...
                permanentMagnet,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');

Pole.compile;
assert(numel(Pole.Regions) == 2);
assert(strcmp(Pole.Regions(1).Name,'pm'));
assert(strcmp(Pole.Regions(2).Name,'DefaultRegion'));
areaDifference = abs(    pi*(0.085^2-0.05^2)/8 ...
                        -Pole.Regions(1).Geometry.area ...
                        -Pole.Regions(2).Geometry.area);
assert(areaDifference < sqrt(eps)*2*pi/80);
assert(area(Pole.Regions(1).Geometry*Pole.Regions(2).Geometry) < sqrt(eps)*2*pi/80);

Pole.removeRegion('pm');
assert(numel(Pole.Regions) == 1);
Pole.compile;

crookedPM = rotate(permanentMagnet,'Rotation',pi/100,'Position',pmPosition);
Pole.addRegion( 'cpm',...
                crookedPM,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');
Pole.compile;
            
assert(Pole.SolutionSpatialFrequency  == 4);
assert(Pole.GeometryFrequency         == 8);
assert(Pole.SolutionHalfWaveSymmetry  == true);
assert(Pole.GeometryMirrorSymmetry    == false);

assert(Pole.ModelSpatialFrequency     == 4);
assert(Pole.ModelHalfWaveSymmetry     == true);

Pole.addRegion( 'pm',...
                permanentMagnet,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');
assert(numel(Pole.Regions) == 3);

Pole.removeRegion( 'cpm');
assert(numel(Pole.Regions) == 2);
            
Parameters.new('fbMinBridgeWidth', 0.0004);

xFluxBarrier1 = pmOuterRadius-pmWidth;
yFluxBarrier1 = -pmLength/2;
bFluxBarrier1 = yFluxBarrier1-1/(tan(pi/nPoles))*xFluxBarrier1;

yFluxBarrier2 = yFluxBarrier1;
xFluxBarrier2 = sqrt((outerRadius-fbMinBridgeWidth)^2-yFluxBarrier2^2);

thetaFluxBarrier3 = asin(fbMinBridgeWidth/outerRadius)-pi/nPoles;
xFluxBarrier3 = (outerRadius-fbMinBridgeWidth)*cos(thetaFluxBarrier3);
yFluxBarrier3 = (outerRadius-fbMinBridgeWidth)*sin(thetaFluxBarrier3);
bFluxBarrier3 = yFluxBarrier3+tan(pi/nPoles)*xFluxBarrier3;

xFluxBarrier4 = (bFluxBarrier3-bFluxBarrier1)/(tan(pi/nPoles)+1/tan(pi/nPoles));
yFluxBarrier4 = -tan(pi/nPoles)*xFluxBarrier4+bFluxBarrier3;

fluxBarrier1Points = [  xFluxBarrier4 yFluxBarrier4;
                        xFluxBarrier3 yFluxBarrier3;
                        xFluxBarrier2 yFluxBarrier2;
                        xFluxBarrier1 yFluxBarrier1]; %polygon points must be
                                                      %specified in a right
                                                      %handed fashion
                    
fluxBarrier1  = Geometry2D.draw(    'Polygon2D',...
                                    'Points',fluxBarrier1Points,...
                                    'PlotStyle',{'w'});
                                
Pole.addRegion( 'fb1',...
                fluxBarrier1,...
                Air,...
                'Static',...
                'Current');
assert(numel(Pole.Regions) == 3);

Pole.compile;
assert(Pole.SolutionSpatialFrequency  == 4);
assert(Pole.GeometryFrequency         == 8);
assert(Pole.SolutionHalfWaveSymmetry  == true);
assert(Pole.GeometryMirrorSymmetry    == false);

assert(Pole.ModelSpatialFrequency     == 4);
assert(Pole.ModelHalfWaveSymmetry     == true);


fluxBarrier2Points = [  xFluxBarrier1 -yFluxBarrier1;
                        xFluxBarrier2 -yFluxBarrier2;
                        xFluxBarrier3 -yFluxBarrier3;
                        xFluxBarrier4 -yFluxBarrier4]; %polygon points must be
                                                       %specified in a right
                                                       %handed fashion

fluxBarrier2  = Geometry2D.draw(    'Polygon2D',...
                                    'Points',fluxBarrier2Points,...
                                    'PlotStyle',{'w'});
                                
Pole.addRegion( 'fb2',...
                fluxBarrier2,...
                Air,...
                'Static',...
                'Current');
assert(numel(Pole.Regions) == 4);

Pole.compile;
assert(Pole.SolutionSpatialFrequency  == 4);
assert(Pole.GeometryFrequency         == 8);
assert(Pole.SolutionHalfWaveSymmetry  == true);
assert(Pole.GeometryMirrorSymmetry    == true);

assert(Pole.ModelSpatialFrequency     == 4);
assert(Pole.ModelHalfWaveSymmetry     == true);


Pole.buildPitch(1/4);
assert(numel(Pole.Regions) == 7);
Pole.buildPitch(1/2);
assert(numel(Pole.Regions) == 13);
Pole.buildPitch(1);
assert(numel(Pole.Regions) == 25);
assert(Pole.GeometryMirrorSymmetry);

magnetIndices = (0:7)*3+1;
magnetAngles  = [0 -3*pi/4 pi/2 -pi/4 pi pi/4 -pi/2 3*pi/4]; 
M  = zeros(8,1);
Mx = zeros(8,1);
My = zeros(8,1);
for i = 1:8
    [M(i),Mx(i),My(i)] = Pole.Regions(magnetIndices(i)).Material.remnantMagnetization;
    ang = atan2(My(i),Mx(i));
    if abs(abs(ang)-pi) < sqrt(eps)
        ang = pi;
    end
    assert(abs(ang-magnetAngles(i)) < sqrt(eps),'%d',i);
end

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