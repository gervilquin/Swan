%% ACADEMIC TEST 0 %%

x0 = [3;2];
j.cF = @(x) x(1)^2 + 2 * x(2)^2;
j.gF = @(x) [2*x(1);4*x(2)];
c.cF = @(x) [2*x(1) + x(2)-9; x(1) + 2*x(2)-10];
c.gF = @(x) [2,1;1,2]';
s.type                       = "IPM";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                     % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND
nConstr                      = 2;                                        % NUMBER OF CONSTRAINTS
s.incrementalScheme.iStep    = 1;
s.incrementalScheme.nSteps   = 1;
s.dualVariable               = [];
s.maxIter                    = [];
s.targetParameters.constr_tol = 1e-3;
s.constraintCase             = {'INEQUALITY','EQUALITY'};
s.maxIter                    = 100;
s.postProcessSettings.shallPrint = false;