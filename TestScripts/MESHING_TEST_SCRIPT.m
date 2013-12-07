clear all;
close all;

display('MESHING_TEST_SCRIPT');

tic
%% Initialize the toolbox
RootObject = MotorProto;

%% Get Parameters
Parameters = RootObject.Parameters;

%% define general pole parameters
Parameters.new('nPoles',            8   );
Parameters.new('outerRadius',       0.085);
Parameters.new('innerRadius',       0.05);

myModel              = RootObject.Model;
Pole                 = myModel.newAssembly('myPole','PoleComponent');
Pole.DefaultMaterial = IronExampleMaterial;
Pole.Poles           = nPoles;
Pole.InnerRadius     = innerRadius;
Pole.OuterRadius     = outerRadius;
Pole.compile;

myMesh = myModel.newMeshFactory;
myMesh.build;
myMesh.plot;