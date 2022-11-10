%% Testing P2 functions
% Create a P1 mesh

clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
femU = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';
mesh = s.mesh;

%% P2 Function
% z.mesh  = mesh;
% z.type    = s.mesh.type;
% z.fValues = femU;
% uP2 = P2Function(z);


%% Projection to P2
% Analytical function
sAF.fHandle = @(x) x(1,:,:).^2;
sAF.ndimf   = 1;
sAF.mesh    = mesh;
uAF = AnalyticalFunction(sAF);

% P1 Function
z.connec  = mesh.connec;
z.type    = mesh.type;
z.fValues = femU;
uP1 = P1Function(z);

% Projector
pp2.mesh   = mesh;
pp2.connec = mesh.connec;
projP2 = Projector_toP2(pp2);
uP2 = projP2.project(uAF);
% uP2.plot(mesh)
% title('P1 (quad linear)')

%% Projection from P2 to P1
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
resp1 = projP1.project(uP2);


%% Stuff to be addressed
% - fValues: how is it passed to P2Function?
% - Do we pass the whole problem mesh? Do we
%