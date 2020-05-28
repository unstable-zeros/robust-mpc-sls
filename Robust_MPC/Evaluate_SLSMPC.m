function [sol, solvertimeTotal] = Evaluate_SLSMPC(sls_mpc, horizon, x0)
%EVALUATE_SLSMPC solve the robust SLS MPC with given horizon and initial
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
range.gamma.lb = 0.01; range.gamma.ub = 10;
range.beta.lb = 0.01; range.beta.ub = 10;
range.tau.lb = 0.01; range.tau.ub = 10;
bisectOptions = struct;
bisectOptions.tol = 1e-2;
[bounds, isFeasible, solvertimeBisect] = sls_mpc.BisectParams(range, bisectOptions);

if isFeasible == 0
    keyboard;
end

%% grid search for feasible hyperparameters
gridDim = struct;
gridDim.num_beta = 5;
gridDim.num_gamma = 5;
gridDim.num_tau = 5;
options = struct;
options.verbose = 0;
[feasibleParams, solvertimeGrid] = sls_mpc.GridSearchParams(bounds, gridDim, options);

params.alpha = 0.5; params.beta = 1.6752; params.gamma = 0.3842; params.tau = 0.2162;

if isempty(feasibleParams)
   keyboard; 
end

%% solve the robust SLS MPC quadratic program with the found hyperparameters
[sol] = sls_mpc.SolveSLSMPC(feasibleParams);
sol.bounds = bounds;
solvertimeMPC = sol.solution.solvertime;

% record the total solvertime 
solvertimeTotal = solvertimeMPC + solvertimeBisect + solvertimeGrid;

sls_mpc.T = original_T;
sls_mpc.x0 = original_x0;

% trajSet = [];
% if options.plot
%     num_traj = options.num_traj;
%     [trajSet, nominalTraj] = sls_mpc.SimulateTraj(num_traj, options);
% end

end

