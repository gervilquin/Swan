filename='RVE_Square_Triangle';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'enforceCh_CCstar_L2','perimeter'};%enforceCh_CCstar_L2
weights=[1 0.05];
constraint = {'volumeConstraint'};
optimizer = 'MMA'; incrementFactor = 1;
filterType = 'P1';

nsteps =10;
Vfrac_final = 0.5;
Perimeter_target=3;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
selectiveC_Cstar = 'Seba';

%  0.0764   -0.0174    0.0088
%    -0.0174    0.0750   -0.0094
%     0.0088   -0.0094    0.0047
% 
% Elapsed time is 58.685115 seconds.