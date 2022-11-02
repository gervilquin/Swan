%% Analytical thermal problem
% syms f(x,y)
% f(x,y) = x^2+y^2;
% fsurf(f)
file = 'Cantileverbeam_Quadrilateral_Bilinear';
a.fileName = file;
s = FemDataContainer(a);

tAF.fHandle = @(x) [ x(1,:,:).^2 + x(2,:,:).^2 ];
tAF.ndimf   = 1;
tAF.mesh    = s.mesh;
xFun = AnalyticalFunction(tAF);

% Quadrature
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');

% Projector to P1
pp1.mesh   = s.mesh;
pp1.connec = s.mesh.connec;
projP1 = Projector_toP1(pp1);
p1fun = projP1.project(xFun);

% Grad
grad = p1fun.computeGradient(quad, s.mesh);

%% 
clc; clear; close all;

% file = 'test2d_triangle';
% file = 'test2d_quad';
% file = 'test3d_hexahedra';
file = 'Cantileverbeam_Quadrilateral_Bilinear';
a.fileName = file;
s = FemDataContainer(a);
fem = FunThermalProblem(s);
fem.solve();

%% Boundary conditions
% logicalfun = @(x) x(1,:,:) == 0;
% logicalfun = @(x) x(1,:,:) >= x(2,:,:);
% logicalfun = @(x) (x(1,:,:)-0.5).^2+(x(2,:,:)-0.5).^2 - 0.3^2 <= 0;
logicalfun = @(x) (x(1,:,:)>=0.1) & (x(1,:,:)<=0.2 & x(2,:,:)>=0.1) & (x(2,:,:)<=0.2);
coor = s.mesh.coord';

z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = double(logicalfun(coor)');
bcMap = P1Function(z);
bcMap.plot(s.mesh)