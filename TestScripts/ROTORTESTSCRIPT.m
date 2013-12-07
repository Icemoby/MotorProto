clear all;close all;
tic;
%% Initialize the toolbox
RootObject = MotorProto('Rotor Model Tutorial');

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general stator parameters
Parameters.new('nPoles',            8   );
Parameters.new('outerRadius',       0.085);
Parameters.new('innerRadius',       0.05);

%% assign rotor parameters;
Rotor              = RootObject.Rotor;
Rotor.Poles        = nPoles;
Rotor.InnerRadius  = innerRadius;
Rotor.OuterRadius  = outerRadius;

%% assign rotor material properties
Rotor.LaminationMaterial = LaminatedIronMaterial;

%% assign mesh parameters
Rotor.MaximumBoundaryEdgeLength = 2*outerRadius*sin(pi/288);
Rotor.InternalEdgeUniformity    = 1.1;

%% define permanent magnet geometry
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
                                        'Material', PermanentMagnetMaterial,...
                                        'PlotStyle',{'m'});

Rotor.add(permanentMagnet);

%% define flux barrier geometry
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
                                    'Material',Air,...
                                    'PlotStyle',{'w'});
                                
Rotor.add(fluxBarrier1);

fluxBarrier2Points = [  xFluxBarrier1 -yFluxBarrier1;
                        xFluxBarrier2 -yFluxBarrier2;
                        xFluxBarrier3 -yFluxBarrier3;
                        xFluxBarrier4 -yFluxBarrier4]; %polygon points must be
                                                       %specified in a right
                                                       %handed fashion

fluxBarrier2  = Geometry2D.draw(    'Polygon2D',...
                                    'Points',fluxBarrier2Points,...
                                    'Material',Air,...
                                    'PlotStyle',{'w'});
                                
Rotor.add(fluxBarrier2);

%% plot
figure;
Rotor.plot;
pause(1);

Rotor.generateMesh;
figure;
Rotor.Mesh.plot
toc