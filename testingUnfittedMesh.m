%% Testing unfitted mesh

clear
clc

file = 'SquareForAniTests';

% Background and boundary mesh
s.isBackgroundMeshRectangularBox = true;
s.inputFile = file;
backBound = BackgroundAndBoundaryMeshCreatorFromInputFile(s);

% Unfitted mesh = mesh
s = [];
s.backgroundMesh = backBound.backgroundMesh;
s.boundaryMesh   = backBound.boundaryMesh;
unfMesh = UnfittedMesh(s);

% LevelSet example
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
sAF.fHandle = @(x) (x(1,:,:)-0.5).^2+(x(2,:,:)-0.5).^2-0.3^2;
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
resAFtoP1 = projP1.project(xFun);

% Unfitted mesh = levelSet
unfMesh.compute(resAFtoP1.fValues);
unfMesh.createPlotter;
% unfMesh.plot;

% Projector with unfitted mesh
s         = [];
s.mesh    = unfMesh.backgroundMesh;
s.fxy     = @(x,y) (x-0.5)^2+(y-0.5)^2-0.3^2;
circleFun = CharacteristicFunction(s);
s = [];
s.mesh   = unfMesh;
s.connec = unfMesh.backgroundMesh.connec;
s.projectorType = 'toH1P1';
projH1P1 = Projector.create(s);
resUnf = projH1P1.project(circleFun);

% Projector with std mesh
s = [];
s.mesh   = unfMesh.backgroundMesh;
s.connec = unfMesh.backgroundMesh.connec;
s.projectorType = 'toH1P1';
projH1P1 = Projector.create(s);
res = projH1P1.project(circleFun);

% 1D - PLOT in function of theta; value(nodes)=f(x)
%    where nodes: 0-eps <= y-mx <= 0+eps
X = unfMesh.backgroundMesh.coord(:,1)-0.5;
Y = unfMesh.backgroundMesh.coord(:,2)-0.5;
theta = 0;
slope = tan(theta);
eps = unfMesh.backgroundMesh.computeMeanCellSize();
fun = Y-slope*X;
k1 = fun >= 0-2*eps;
k2 = fun <= 0+2*eps;
nodesk = find(k1==k2);
l = sqrt(X(nodesk).^2+Y(nodesk).^2);
figure
plot(l,resUnf.fValues(nodesk),'.')
hold on
plot(l,res.fValues(nodesk),'.')
hold off
grid on
grid minor
xlabel('Length','Interpreter','latex')
ylabel('CharFun','Interpreter','latex')
legend('Unfitted mesh','Standard mesh','Interpreter','latex')
ylim([-0.1 1.1])

figure
semilogy(l,abs(resUnf.fValues(nodesk)-res.fValues(nodesk)),'.')
grid on
grid minor
xlabel('Length','Interpreter','latex')
ylabel('Absolute error','Interpreter','latex')

% Volume computation
Vreal = pi*0.3^2;

quad = Quadrature.set('TRIANGLE');
quad.computeQuadrature('LINEAR');

VunfInner = unfMesh.innerMesh.mesh.computeDvolume(quad);
VunfInnerCut = unfMesh.innerCutMesh.mesh.computeDvolume(quad);
Vunf = sum(VunfInner)+sum(VunfInnerCut);