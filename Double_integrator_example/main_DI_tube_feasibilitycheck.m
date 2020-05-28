%% The following codes generate Fig.4 in https://arxiv.org/pdf/1911.06842.pdf
clear;

%% construct robust MPC
load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

sls_mpc = SLSDesign( sysdata );

%% find the disturbance invariant set Z_inv
options = struct;
uncertain_sys = UncertainSystem(sysdata, options);
[Z_inv, isConverge] = uncertain_sys.minInvSet(10);

% construction of unit inf-norm ball Z_unit
Zvertices = [-1 -1; -1 1; 1 -1; 1 1];
Z_unit = Polyhedron(Zvertices);
Z_unit.H; % initialize the H representation of Z_unit

%% initial condition feasibility evaluation of tube MPC
% extract the robust maximum invariant set RIS
RIS = sls_mpc.terminalConstraints;

%  grid search inside RIS
samples = RIS.grid(13);
x0Set = samples;

tubeOpts = struct;
tubeOpts.verbose = 2;
tubeOpts.normType = 2;

% verify feasibility of tube MPC on the given set of initial
% conditions. Change the tube shape parameter to either Z_inv or Z_unit
[solCell, classSet] = VerifyFeasibility_Tube(sls_mpc, x0Set, Z_inv, tubeOpts);

figure;
RIS.plot();
hold on
fset = classSet.feasibleSet; % feasible x0
ifset = classSet.infeasibleSet; % infeasible x0
scatter(fset(:,1), fset(:,2), 'ys');
scatter(ifset(:,1), ifset(:,2), 'go');
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title('Feasibility check of tube MPC', 'Interpreter', 'Latex', 'FontSize', 20);
legend('max RIS', 'feasible', 'infeasible');

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

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);