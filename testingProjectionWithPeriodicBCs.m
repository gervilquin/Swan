% TESTING PROJECTION WITH PERIODIC BCs
%% Get mesh
clear; close all;

file = 'RVE_Square_Triangle';
a.fileName = file;
m = FemDataContainer(a);

%% Create characteristic function
s.mesh    = m.mesh;
s.fxy     = @(x,y) ((x-0.5)^2+y^2-0.2^2).*(y<=0.2) + ...
                   (x^2+(y-0.7)^2-0.2^2).*(x<=0.2) + ...
                   ((x-0.4)^2+(y-1)^2-0.2^2).*(y>=0.8) + ...
                   ((x-1)^2+(y-0.4)^2-0.2^2).*(x>=0.8) + 0.001;
circleFun = CharacteristicFunction(s);

%% Projection

% CharFun to P1 Function in H1 domain
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
x.projectorType = 'toPeriodicH1P1';
x.masterSlave = m.mesh.masterSlaveNodes;
projH1P1 = Projector.create(x);
resCharFunInH1P1 = projH1P1.project(circleFun);
resCharFunInH1P1.plot(m.mesh);