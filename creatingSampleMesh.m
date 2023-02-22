clc; clear;

N = 10;
coordX = linspace(0,1,N);
coordY = linspace(0,1,N);
[X,Y] = meshgrid(coordX,coordY);
Z     = zeros(size(X));
[conn, coordV] = mesh2tri(X,Y,Z, 'f');

s.coord  = coordV(:,1:2);
s.connec = conn;
mesh = Mesh(s);
mesh.computeMasterSlaveNodes();

nNodes = size(s.coord, 1);
coordsDef = zeros(nNodes, 4);

for i = 1:size(s.coord, 1)
    coordsDef(i, 1) = i;
end

coordsDef(:, 2:3) = s.coord;

nElem = size(conn, 1);

connDef = zeros(nElem, 4);

connDef(:, 2:4) = conn;
for i=1:nElem
    connDef(i, 1) = i;
end
