%% generate Fig.3 in https://arxiv.org/pdf/1911.06842.pdf
% compare the solver time of SLS MPC and tube MPC by running on a robust
% MPC exmaple over a range of horizons.
clear;

%% record the SLS MPC solver time
load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

sls_mpc = SLSDesign( sysdata );

horizonList = [4 8 12 16 20];
num_horizon = length(horizonList);

sls_mpc_timeList = zeros(1, num_horizon);
solList = cell(1, num_horizon);
for ii = 1:num_horizon
    fprintf('Running time test: %d out of %d\n', ii, num_horizon);
    horizon = horizonList(ii);
    [sol, solvertimeTotal] = Evaluate_SLSMPC(sls_mpc, horizon);
    solList{ii} = sol;
    sls_mpc_timeList(ii) = solvertimeTotal;
end
figure; 
plot(horizonList, sls_mpc_timeList,'bo-');
xlabel('horizon', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('solver time', 'Interpreter', 'Latex', 'FontSize', 18);
title('running time of SLS MPC', 'Interpreter', 'Latex', 'FontSize', 20);

%% record the tube MPC solver time
% construct tube mpc model
load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

sls_mpc = SLSDesign( sysdata );

% compute minimum invariant set
options = struct;
uncertain_sys = UncertainSystem(sysdata, options);
optRIS.robust = 1;
optRIS.minVol = 0.5;
[Z_inv, isConverge] = uncertain_sys.minInvSet(10);

Z = Z_inv;

tubeOpts = struct;
tubeOpts.verbose = 2;
tubeOpts.normType = 2;

% compute the running time
num_horizon = length(horizonList);
tubempc_timeList = zeros(1, num_horizon);
solList = cell(1, num_horizon);
for ii = 1:num_horizon
    fprintf('%d out of %d\n', ii, num_horizon);
    sls_mpc.T = horizonList(ii);
    [sol] = sls_mpc.SolveTubeMPC(Z, tubeOpts);
    solList{ii} = sol;
    tubempc_timeList(ii) = sol.diagnostics.solvertime;
end
figure; 
plot(horizonList, tubempc_timeList, 'bo-');
xlabel('horizon', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('cpu time', 'Interpreter', 'Latex', 'FontSize', 18);
title('running time of tube mpc', 'Interpreter', 'Latex', 'FontSize', 20);

%% compare the solver time of SLS MPC and tube MPC
figure;
plot(horizonList, tubempc_timeList,'bo-', 'LineWidth', 2);
hold on
plot(horizonList, sls_mpc_timeList,'rs-', 'LineWidth', 2);

xlabel('horizon', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('solver time/seconds', 'Interpreter', 'Latex', 'FontSize', 18);
title('running time of tube/sls mpc', 'Interpreter', 'Latex', 'FontSize', 20);
legend('tube', 'sls');

