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