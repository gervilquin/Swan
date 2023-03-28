%% Plot results

% Plot mesh with some nodes
nodesL = obj.bc.masterSlave(:,1);
coordsL = obj.mesh.coord(nodesL, :);

obj.mesh.plot;
plot(coordsL(:,1), coordsL(:,2), '+');

% P1Functions
a.mesh = obj.mesh;
fvrshp = reshape(u, [2 265]);
a.fValues = zeros(obj.mesh.nnodes,1);
a.fValues = fvrshp';

p1f = P1Function(a);
p1f.plot;
