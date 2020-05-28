function [sol] = Evaluate_CoarseSLSMPC(sls_mpc, horizon, x0)
%EVALUATE_COARSESLSMPC solve the robust coarse SLS MPC with given horizon and initial
%condition following the 'bound hyperparamters -> grid search feasible
%hyperparamters -> solve QP' pipeline.

if nargin < 2
    horizon = sls_mpc.T;
    x0 = sls_mpc.x0;
end

if nargin < 3
    x0 = sls_mpc.x0;
end

original_T = sls_mpc.T;
original_x0 = sls_mpc.x0;

sls_mpc.T = horizon;
sls_mpc.x0 = x0;

%% find hyperparameter bounds
range = struct;
range.tau.lb = 0.01; range.tau.ub = 10;
bisectOptions = struct;
bisectOptions.tol = 1e-2;
[bounds, isFeasible] = sls_mpc.CoarseBisectParams(range, bisectOptions);

if isFeasible == 0
    keyboard;
end

%% grid search for feasible hyperparameters
gridDim = struct;
gridDim.num_tau = 3;
options = struct;
options.verbose = 0;
[feasibleParams] = sls_mpc.CoarseGridSearchParams(bounds, gridDim, options);

params.alpha = 0.5; params.beta = 1.6752; params.gamma = 0.3842; params.tau = 0.2162;

if isempty(feasibleParams)
   keyboard; 
end

%% solve the robust coarse SLS MPC quadratic program with the found hyperparameters
[sol] = sls_mpc.SolveCoarseSLSMPC(feasibleParams);
sol.bounds = bounds;

sls_mpc.T = original_T;
sls_mpc.x0 = original_x0;



end

