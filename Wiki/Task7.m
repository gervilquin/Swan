%% Swan - Wiki
% The student's guide to clean code development
% Task 7: UML of FEM code in Swan repository

% Instructions: run the following code, selecting previously the 'Swan'
% main folder as your current matlab path

file = 'test2d_triangle_non_null';
a.fileName = file;
s = FemDataContainer(a);
s = ObjectSaver.saveObj(s);
s.builderType = 'MONOLITIC';
femM = FEM.create(s);
femM.solve();

file = 'test2d_triangle_non_null';
a.fileName = file;
s = FemDataContainer(a);
s = ObjectSaver.saveObj(s);
s.builderType = 'REDUCED';
femR = FEM.create(s);
femR.solve();

err = femM.variables.d_u- femR.variables.d_u;
norm(err(:))
