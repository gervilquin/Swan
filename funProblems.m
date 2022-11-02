%% Testing funcitons for FEM
clc; clear; close all;

file = 'test2d_triangle';
% file = 'test2d_quad';
% file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);

% Boundary conditions
nP.type      = 'Neumann';
nP.value     = [0, -1]; % that zero... meh
nP.domainFun = @(x) (x(1,:,:) == 0.04) & (x(2,:,:) == 0.02);
pointF = BoundaryCondition(nP);

dP.type      = 'Dirichlet';
dP.value     = [0, 0];
dP.domainFun = @(x) x(1,:,:) == 0;
fixedU = BoundaryCondition(dP);

s.bc.dirichlet = fixedU;
s.bc.neumann   = pointF;

% Anyway
fem = FunElasticProblem(s);
fem.solve();

%% Boundary conditions as functions
% AnalyticalFunction
sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
% sAF.fHandle = @(x) [cos(x(1,:,:).*x(2,:,:)); x(1,:,:).*x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = s.mesh;
xFun = AnalyticalFunction(sAF);