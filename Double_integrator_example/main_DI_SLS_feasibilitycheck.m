clear;
%% construct robust SLS MPC
load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

% initialize the robust SLS MPC problem
sls_mpc = SLSDesign( sysdata );

x0 = sls_mpc.x0;

% solve the robust SLS MPC at the initial condition x0
[sol, solvertimeTotal] = Evaluate_SLSMPC(sls_mpc);
% extract the lower and upper bounds on the hyperparameters with given x0
% since the lower bound on tau and beta are independent of x0
bounds = sol.bounds;

%% initial condition feasibility evaluation 
% we do a grid search of the initial conditions over the maximum robust
% invariant set and test feasibility of the robust SLS MPC formulation. The
% following codes generate Fig.5(a) in https://arxiv.org/pdf/1911.06842.pdf

% extract the maximum robust invariant set RIS
RIS = sls_mpc.terminalConstraints;

% grid search inside RIS
samples = RIS.grid(13);
x0Set = samples;

feasibility_opt = struct;

gridDim = struct;
% grid search resolution 3 by 3 by 3
gridDim.num_beta = 3;
gridDim.num_gamma = 3;
gridDim.num_tau = 3;

gridparams_opt = struct;
gridparams_opt.verbose = 0;

% grid search range
range = struct;
range.gamma.lb = 0.01; range.gamma.ub = 10;
range.beta.lb = 0.01; range.beta.ub = 10;
range.tau.lb = 0.01; range.tau.ub = 10;

bisect_opt = struct;
bisect_opt.init.tau.lb = bounds.tau.lb;
bisect_opt.init.beta.lb = bounds.beta.lb;
bisect_opt.tol = 1e-2;
bisect_opt.verbose = 0;

feasibility_opt.range = range;
feasibility_opt.gridparams_opt = gridparams_opt;
feasibility_opt.gridDim = gridDim;
feasibility_opt.bisect_opt = bisect_opt;

% verify feasibility of robust SLS MPC on the given set of initial
% conditions
[solCell, sets] = VerifyFeasibility(sls_mpc, x0Set, feasibility_opt);

fset = sets.feasibleSet; % feasible x0
ifset = sets.infeasibleSet; % infeasible x0
uvset = sets.unverifiedSet; % unverified x0

%% plot the initial conditions based on feasibility
figure;
color_1 = [178 102 255]./255;
PolytopePlot.show_convex(sls_mpc.stateConstraints, color_1, 'LineStyle', 'none');
hold on
color_2 = [0 153 76]./255;
if ~isempty(sls_mpc.terminalConstraints)
    PolytopePlot.show_convex(sls_mpc.terminalConstraints, color_2, 'LineStyle', 'none');
end
scatter(fset(:,1), fset(:,2), 'ys', 'filled');
scatter(ifset(:,1), ifset(:,2), 'ro', 'filled');
scatter(uvset(:,1), uvset(:,2), 'w*');
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);