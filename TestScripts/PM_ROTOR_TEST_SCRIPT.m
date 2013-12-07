clear all;
close all;

display('PM_ROTOR_TEST_SCRIPT');
tic
%% Initialize the toolbox
RootObject = MotorProto('Permanent Magnet Sychronous Rotor Demo');

%% Test adding a StatorLamination object to the Model
Model = RootObject.Model;
pole  = Model.newComponent('myPole','PoleComponent');

assert(pole == Model.Components(1));
assert(strcmp(class(Model.Components(1)),'PoleComponent'));
assert(strcmp(Model.Components(1).Name,'myPole'));

%% Test removing a StatorLamination methods
Model.removeComponent('myPole');
assert(isempty(Model.Components));
assert(strcmp(pole.Name,'myPole'));

%% Test removing a StatorLamination which is not there
Model.removeComponent('myPole');
assert(isempty(Model.Components));

%% Try adding two components with the same name
P1 = Model.newComponent('myPole','PoleComponent');
assert(numel(Model.Components)==1);
assert(strcmp(Model.Components(1).Name,'myPole'));

try
    Model.newComponent('myPole','PoleComponent');
end
assert(numel(Model.Components)==1);
assert(strcmp(Model.Components(1).Name,'myPole'));
assert(Model.Components(1) == P1);

%% Try adding a component with a different name
P2 = Model.newComponent('myOtherPole','PoleComponent');
assert(numel(Model.Components)==2);
assert(strcmp(Model.Components(2).Name,'myOtherPole'));

%% Try removing a component that's not there
Model.removeComponent('myPole_1');
assert(isvalid(P1));
assert(isvalid(P2));
assert(numel(Model.Components) == 2);
assert(strcmp(Model.Components(1).Name,'myPole'));

%% Try deleting components
Model.deleteComponent('myOtherPole');
assert(~isvalid(P2));
assert(numel(Model.Components) == 1);

%% Try creating a new component with a non-string name
try
	P3 = Model.newComponent(1,'PoleComponent');
catch ME
    P3 = [];
end
assert(isempty(P3));
assert(numel(Model.Components) == 1);
assert(strcmp(Model.Components(1).Name,'myPole'));

P2 = Model.newComponent('1','PoleComponent');
assert(numel(Model.Components) == 2);
assert(strcmp(Model.Components(2).Name,'1'));

%% Try adding a component which has a null name
P2 = PoleComponent;
try 
    Model.addComponent(P2);
end
assert(numel(Model.Components) == 2);
Model.deleteComponent('1');
assert(numel(Model.Components) == 1);

pole = Model.Components(1);

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general pole parameters
Parameters.new('nPoles',            8   );
Parameters.new('outerRadius',       0.085);
Parameters.new('innerRadius',       0.05);

pole.Poles        = nPoles;
pole.InnerRadius  = innerRadius;
pole.OuterRadius  = outerRadius;

assert(pole.SolutionSpatialFrequency  == 4);
assert(pole.GeometryFrequency         == 8);
assert(pole.SolutionHalfWaveSymmetry  == true);
assert(pole.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency == 4);
assert(Model.HasHalfWaveSymmetry      == true);

pole.build;
assert(pole.SolutionSpatialFrequency  == 4);
assert(pole.GeometryFrequency         == 8);
assert(pole.SolutionHalfWaveSymmetry  == true);
assert(pole.GeometryMirrorSymmetry    == true);
assert(Model.SolutionSpatialFrequency == 4);
assert(Model.HasHalfWaveSymmetry      == true);

%% Try assigning some pole properties
pole.DefaultMaterial = IronExampleMaterial;
try
    pole.DefaultMaterial = 1;
end
assert(isa(pole.DefaultMaterial,'IronExampleMaterial'));

pole.build;
assert(numel(pole.Regions) == 1);
assert(strcmp(pole.Regions.Name,'myPole_DefaultRegion'));
assert(isa(pole.Regions.Geometry,'Sector'));
assert(all(pole.Regions.Geometry.Radius.Value == [0.05 0.085]));
assert(abs(pole.Regions.Geometry.Angle.Value - 2*pi/8) < eps*2*pi);
assert(pole.Regions.Geometry.Rotation == -pi/8);
assert(isa(pole.Regions.Material,'IronExampleMaterial'));

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

pole.addRegion( 'pm',...
                permanentMagnet,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');

pole.build;
assert(numel(pole.Regions) == 2);
assert(numel(pole.Components) == 2);
assert(strcmp(pole.Regions(1).Name,'pm'));
assert(strcmp(pole.Regions(2).Name,'myPole_DefaultRegion'));
assert(strcmp(pole.Components(1).Name,'pm'));
assert(strcmp(pole.Components(2).Name,'myPole_DefaultRegion'));
areaDifference = abs(    pi*(0.085^2-0.05^2)/8 ...
                        -pole.Regions(1).Geometry.area ...
                        -pole.Regions(2).Geometry.area);
assert(areaDifference < sqrt(eps)*2*pi/80);
assert(area(pole.Regions(1).Geometry*pole.Regions(2).Geometry) < sqrt(eps)*2*pi/80);

pole.removeComponent('pm');
assert(numel(pole.Regions) == 1);
assert(numel(pole.Components) == 1);
pole.build;

crookedPM = rotate(copy(permanentMagnet),pi/100,pmPosition);
pole.addRegion( 'cpm',...
                crookedPM,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');
pole.build;
            
assert(pole.SolutionSpatialFrequency  == 4);
assert(pole.GeometryFrequency         == 8);
assert(pole.SolutionHalfWaveSymmetry  == true);
assert(pole.GeometryMirrorSymmetry    == false);

assert(Model.SolutionSpatialFrequency == 4);
assert(Model.HasHalfWaveSymmetry      == true);

pole.addRegion( 'pm',...
                permanentMagnet,...
                PermanentMagnetExampleMaterial,...
                'Dynamic',...
                'Current');
assert(numel(pole.Regions) == 3);
assert(numel(pole.Components) == 3);

pole.removeComponent( 'cpm');
assert(numel(pole.Regions) == 2);
assert(numel(pole.Components) == 2);
            
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
                                
pole.addRegion( 'fb1',...
                fluxBarrier1,...
                Air,...
                'Static',...
                'Current');
assert(numel(pole.Regions) == 3);
assert(numel(pole.Components) == 3);

pole.build;
assert(pole.SolutionSpatialFrequency  == 4);
assert(pole.GeometryFrequency         == 8);
assert(pole.SolutionHalfWaveSymmetry  == true);
assert(pole.GeometryMirrorSymmetry    == false);

assert(Model.SolutionSpatialFrequency == 4);
assert(Model.HasHalfWaveSymmetry      == true);


fluxBarrier2Points = [  xFluxBarrier1 -yFluxBarrier1;
                        xFluxBarrier2 -yFluxBarrier2;
                        xFluxBarrier3 -yFluxBarrier3;
                        xFluxBarrier4 -yFluxBarrier4]; %polygon points must be
                                                       %specified in a right
                                                       %handed fashion

fluxBarrier2  = Geometry2D.draw(    'Polygon2D',...
                                    'Points',fluxBarrier2Points,...
                                    'PlotStyle',{'w'});
                                
pole.addRegion( 'fb2',...
                fluxBarrier2,...
                Air,...
                'Static',...
                'Current');
assert(numel(pole.Regions) == 4);
assert(numel(pole.Components) == 4);

pole.build;
assert(pole.SolutionSpatialFrequency  == 4);
assert(pole.GeometryFrequency         == 8);
assert(pole.SolutionHalfWaveSymmetry  == true);
assert(pole.GeometryMirrorSymmetry    == true);

assert(Model.SolutionSpatialFrequency == 4);
assert(Model.HasHalfWaveSymmetry      == true);


pole.build(1/4);
assert(numel(pole.Regions) == 7);
assert(numel(pole.Components) == 7);
pole.build(1/2);
assert(numel(pole.Regions) == 13);
assert(numel(pole.Components) == 13);
pole.build(1);
assert(numel(pole.Regions) == 25);
assert(numel(pole.Components) == 25);
assert(pole.GeometryMirrorSymmetry);

magnetIndices = (0:7)*3+1;
magnetAngles  = [0 -3*pi/4 pi/2 -pi/4 pi pi/4 -pi/2 3*pi/4]; 
M  = zeros(8,1);
Mx = zeros(8,1);
My = zeros(8,1);

for i = 1:8
    [M(i),Mx(i),My(i)] = pole.Regions(magnetIndices(i)).Material.remnantMagnetization;
    ang(i) = atan2(My(i),Mx(i));
    if abs(abs(ang(i))-pi) < sqrt(eps)
        ang(i) = pi;
    end
end
assert(all(abs(ang-magnetAngles) < sqrt(eps)));

%% Generate the mesh factory
poleMesh = Model.Mesh(1);
poleMesh.MaximumElementSize = (outerRadius - innerRadius) / 20;
poleMesh.MinimumElementSize = 0;
poleMesh.UseUniformGrid     = false;
poleMesh.build;

assert(numel(poleMesh.Boundaries) == 98);
assert(isempty(setxor(unique(poleMesh.ElementRegions),1:25)));

RootObject.configureAlgorithm('Static');

%% test run
s = RootObject.run;

%% test various matrix properties
K = s.Matrices.K(0);
isSymmetric        = max(max(abs(K-K.'))) / max(max(abs(K))) < eps;
[~,p]              = chol(K);
isPositiveDefinite = (p == 0);
assert(isPositiveDefinite);
assert(isSymmetric);

G            = s.Matrices.G(0,s.Results{1});
isSymmetric  = max(max(abs(G-G.'))) / max(max(abs(G))) < eps;
assert(isSymmetric);

J = K + G;
isSymmetric        = max(max(abs(J-J.'))) / max(max(abs(J))) < eps;
[~,p]              = chol(J);
isPositiveDefinite = (p == 0);
assert(isPositiveDefinite);
assert(isSymmetric);
toc 