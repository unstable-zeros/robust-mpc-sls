%% generate Fig.2 in https://arxiv.org/pdf/1911.06842.pdf
% plot trajectories of a robust MPC example through both SLS MPC and tube MPC.
clear;

%% construct MPC
load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

sls_mpc = SLSDesign( sysdata );

horizon = sls_mpc.T;
horizon

x0 = sls_mpc.x0;

%% find SLS MPC solutions
% find robust SLS MPC solution
[sol, solvertimeTotal] = Evaluate_SLSMPC(sls_mpc);

feasibleParams = struct;
feasibleParams.gamma = sol.gamma;
feasibleParams.beta = sol.beta;
feasibleParams.tau = sol.tau;
feasibleParams.alpha = sol.alpha;

% simulate SLS MPC trajectories
num_traj = 10;

options = struct;
options.ismemory = 0;
options.num_traj = 10;
options.plot = 1;

% simulate several trajectories under model uncertainty
[sls_traj_set, nominalTraj] = sls_mpc.SimulateTraj(num_traj, options);

% plot robust SLS MPC trajectories
PlotTraj(sls_mpc, sls_traj_set);

%% simulate tube MPC trajectories

% find disturbance invariant set Z_inv
optUncertain = struct;
uncertain_sys = UncertainSystem(sysdata, optUncertain);
optRIS.robust = 1;
optRIS.minVol = 0.5;
[Z_inv, isConverge] = uncertain_sys.minInvSet(10);
Z = Z_inv;

% find unit inf-norm ball Z_unit
Zvertices = [-1 -1; -1 1; 1 -1; 1 1];
Z_unit = Polyhedron(Zvertices);
Z_unit.H;

tubeOpts = struct;
tubeOpts.verbose = 2;
tubeOpts.normType = 2;

% solve for the tube MPC control law
% set Z = Z_inv or Z = Z_unit
[tube_sol] = sls_mpc.SolveTubeMPC(Z, tubeOpts);

% implement the tube MPC control law on the previous uncertain system
% examples stored in sls_traj_set for comparison
[tubeTrajSet] = TubePlotTraj(sls_mpc, tube_sol, sls_traj_set );
 