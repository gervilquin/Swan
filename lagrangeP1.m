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
figure
obj.mesh.plot;
xlabel('x');
ylabel('y');
% plot(coordsL(:,1), coordsL(:,2), '+')
plot(obj.mesh.coord([2, 3, 207, 208],1), obj.mesh.coord([2, 3, 207, 208],2), '+')


z.coord = coordsL;
z.connec(:,1) = 1:size(coordsL,1)-1;
z.connec(:,2) = 2:size(coordsL,1);
z.kFace = -1;
m = Mesh(z);

quad = Quadrature.set(m.type);
quad.computeQuadrature('LINEAR');
aa.mesh = m;
aa.fValues = L(1:13);
auxV = zeros(13+2, 1);
auxV(2:14) = L(1:13);
auxV(1) = LDir(1);
auxV(15) = LDir(5);
aa.fValues = auxV;
linep1f = P1Function(aa);
valors = linep1f.evaluate(quad.posgp);
dV = m.computeDvolume(quad);
sigma = sum(sum(squeeze(valors).*squeeze(dV)));

%% Plot fluctuations/displacements over the whole mesh
fVrshp = reshape(u, [2 265]);
a.mesh = obj.mesh;
a.fValues = fVrshp';
uFun = P1Function(a);
uFun.plot()

%% Plot lagrange multipliers
% Periodic L plot
% nodesL = obj.bc.masterSlave(1:13,1);
% coordsL = obj.mesh.coord(nodesL,:);
% fV = zeros(obj.mesh.nnodes, 1);
% fV(nodesL) = L(1:13);
% a.mesh = obj.mesh;
% a.fValues = fV;
% lFun = P1Function(a);

fV = zeros(obj.mesh.nnodes*2,1);
dofsM1 = obj.bc.periodic_free(1:13);
dofsM2 = obj.bc.periodic_free(end-12:end);
% fV(dofsM1) = L(1:13);
% fV(dofsM2) = L(end-12:end);
fV(dofsM1) = L(14:26);
% fV(dofsM2) = L(27:39);
fVlambdasRshp = reshape(fV, [2 265]);
a.fValues = fVlambdasRshp';
a.mesh = obj.mesh;
lFun = P1Function(a);

% Dirichlet L plot - OK
dofsD = obj.bc.dirichlet;
sDir = size(obj.bc.dirichlet, 1);
lambdaD = zeros(sDir, 1);
if obj.vstrain(1) == 1
    lambdaD(1) = LDir(1) + LDir(3);
    lambdaD(2) = LDir(4);
    lambdaD(4) = LDir(7);
    lambdaD(5) = LDir(2) + LDir(5);
    lambdaD(6) = LDir(6);
    lambdaD(8) = LDir(8);
elseif obj.vstrain(2) == 1
    lambdaD(1) = LDir(3);
    lambdaD(2) = LDir(1) + LDir(4);
    lambdaD(3) = LDir(5);
    lambdaD(4) = LDir(2) + LDir(6);
    lambdaD(5) = LDir(7);
    lambdaD(7) = LDir(8);
else
    lambdaD(1) = LDir(1);
    lambdaD(5) = LDir(2);
end
% lambdaD = LDir;
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

% Print GiD results: Lagrange multipliers
a.filename = 'LagrangeGiDFirst';
allLp1f.print(a);

% Print GiD results: Fluctuations
b.filename = 'DisplacementsGiDFirst';
fVrshp = reshape(u, [2 265]);
b.mesh = obj.mesh;
b.fValues = fVrshp';
gidPrint = P1Function(b);
gidPrint.print(b);


% PRINT MACRO RESULTS 
b.filename = 'DisplacementsGiDFirst';
fVrshp = reshape(u, [2 13]);
b.mesh = obj.mesh;
b.fValues = fVrshp';
gidPrint = P1Function(b);
gidPrint.print(b);




