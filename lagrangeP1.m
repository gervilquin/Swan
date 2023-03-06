%% How to run this
% - Add a breakpoint at MicroBuilder.solveSystem
% - Run BoundaryCondTestsSuite
% - Run this file
% - Some sections may not run at first, just copy and paste the code at the
%   command window and it *should* work

%% Plot Lagrange multipliers at the left
% also sigma
nodesL = obj.bc.masterSlave(:,1);
coordsL = obj.mesh.coord(nodesL,:);

nodesL = obj.bc.masterSlave(1:13,1);
coordsL = obj.mesh.coord(nodesL,:);
plot(coordsL(:,1), coordsL(:,2), '+')

z.coord = coordsL;
z.connec(:,1) = 1:size(coordsL,1)-1;
z.connec(:,2) = 2:size(coordsL,1);
z.kFace = -1;
m = Mesh(z);

quad = Quadrature.set(m.type);
quad.computeQuadrature('LINEAR');
aa.mesh = m;
aa.fValues = L(1:13);
linep1f = P1Function(aa);
valors = linep1f.evaluate(quad.posgp);
dV = m.computeDvolume(quad);
sigma = sum(sum(squeeze(valors).*squeeze(dV)));

%% Plot fluctuations over the whole mesh
fVrshp = reshape(u, [2 265]);
a.mesh = obj.mesh;
a.fValues = fVrshp';
uFun = P1Function(a);
uFun.plot()

%% Plot lagrangians noch einmal
nodesL = obj.bc.masterSlave(1:13,1);
coordsL = obj.mesh.coord(nodesL,:);
fV = zeros(obj.mesh.nnodes, 1);
fV(nodesL) = L(1:13);
a.mesh = obj.mesh;
a.fValues = fV;
lFun = P1Function(a);

dofsD = obj.bc.dirichlet;
lambdaD = sol(end-7:end);
fVlambdas = zeros(obj.mesh.nnodes*2, 1);
fVlambdas(dofsD) = lambdaD;
fVlambdasRshp = reshape(fVlambdas, [2 265]);
a.mesh = obj.mesh;
a.fValues = fVlambdasRshp';
lambdasp1f = P1Function(a);
lambdasp1f.plot
ltot = lFun.fValues + lambdasp1f.fValues;
clc
a.mesh = obj.mesh;
a.fValues = ltot;
allLp1f = P1Function(a);
allLp1f.plot