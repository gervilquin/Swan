%% Micro
% Get FEM results
clear; close all;

file = 'test2d_micro';
file = 'micro_hole_example';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.computeChomog();
fem.print(file)

%% Macro
% Get FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
fem.print(file)