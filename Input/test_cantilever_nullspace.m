% filename = 'CantileverVertical';
filename = 'MBB';

ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';

initial_case = 'full'; % full OR circleInclusion

cost = {'compliance','anisotropicPerimeter2D'};
weights = [1,1e-5];

constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};

optimizerUnconstrained = 'SLERP'; 
optimizer = 'NullSpace';
incrementFactor = 2;
designVariable = 'LevelSet';
filterType = 'P1';
fracRadius = 0.3;

nsteps = 1;
Vfrac_final = 0.2;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

% For all tests
plotting = false;
printing = true;
printing_physics = false;
monitoring = true;
monitoring_interval = 5;
maxiter = 300; % 2200