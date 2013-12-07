clear all;
close all;
clc

GEOMETRY_TEST_SCRIPT
PM_ROTOR_TEST_SCRIPT
IM_ROTOR_TEST_SCRIPT
STATOR_TEST_SCRIPT
%MESHING_TEST_SCRIPT

%% TODO
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
%   Get rid of the local axis in the Geometry class (tends to be slow). 
%   Creating a "rotatable" mixin might be a good solution. This could replace
%   the local axis object in the material property subclasses as well.
%
%   Vectorize, vectorize, vectorize, minimize the number of function calls
%       -The Geometry1D and Geometry2D "inOn" methods could be vectorized