function [tubeTrajSet] = TubePlotTraj(sls_mpc, tube_sol, sls_traj_set )
%TUBEPLOTTRAJ Plot the tubes and trajectories of the tube MPC. For each
%tube MPC trajectory, the same disturbance and model uncertainty sequences
%are used as in the reference robust SLS MPC exmaple given by sls_traj_set.

num_traj = length(sls_traj_set);

% recover model 
T = sls_mpc.T; x0 = sls_mpc.x0;
n = sls_mpc.stateDim; p = sls_mpc.inputDim;
sigmaW = sls_mpc.sigmaW;

Ahat = sls_mpc.Ahat; Bhat = sls_mpc.Bhat;
epsA = sls_mpc.epsA; epsB = sls_mpc.epsB;

% extract Delta_A , Delta_B and w sequences from sls_traj_set
DeltaA_set = cell(1,num_traj);
DeltaB_set = cell(1,num_traj);
w_set = cell(1, num_traj);
for ii = 1:num_traj
   DeltaA_seq = cell(1, T);
   DeltaB_seq = cell(1, T);
   for jj = 1:T
       DeltaA_seq{jj} = sls_traj_set{ii}.DeltaA_seq{jj}(:, end-n+1:end);
       DeltaB_seq{jj} = sls_traj_set{ii}.DeltaB_seq{jj}(:, end-p+1:end);
   end
   DeltaA_set{ii} = DeltaA_seq;
   DeltaB_set{ii} = DeltaB_seq; 
   w_set{ii} = sls_traj_set{ii}.w_seq;
end
u0 = tube_sol.u0;
tube_seq = tube_sol.tubes;
u_tubes = tube_sol.u_tubes;
vertices_set = tube_sol.vertices_set;

% compute the tube MPC trajectories
tubeTrajSet = cell(1, num_traj);
for ii = 1:num_traj
   traj = struct;
   x_seq = [x0]; u_seq = [u0]; 
   
   DeltaA_seq = DeltaA_set{ii}; DeltaB_seq = DeltaB_set{ii};
   w_seq = w_set{ii};
   
   x_tube = x0; 
   for jj = 1:T
       u_cur = u_seq(:, end);
       dist_A = DeltaA_seq{jj}*x_tube;
       dist_B = DeltaB_seq{jj}*u_cur;
       w = w_seq(:,jj + 1); % skip the initial condition term
         
       x_tube = Ahat*x_tube + dist_A + Bhat*u_cur + dist_B + w;
       
       tube = vertices_set{jj};

       [coef, isfeasible] = FindHullCoeff(x_tube, cell2mat(tube)', 1);
       if isfeasible ~=0
           keyboard;
       end
       % find tube MPC control inputs
       u_next = cell2mat(u_tubes{jj})*coef;
       u_seq = [u_seq u_next];
       x_seq = [x_seq x_tube];
       
       end
   traj.x_seq = x_seq;
   traj.u_seq = u_seq;
   traj.w_seq = w_seq;
   traj.DeltaA_seq = DeltaA_seq;
   traj.DeltaB_seq = DeltaB_seq; 

   tubeTrajSet{ii} = traj;
end

%% plot the tube MPC trajectories
figure;
color_1 = [178 102 255]./255;
PolytopePlot.show_convex(sls_mpc.stateConstraints, color_1, 'LineStyle', 'none');
hold on
color_2 = [0 153 76]./255;
if ~isempty(sls_mpc.terminalConstraints)
    PolytopePlot.show_convex(sls_mpc.terminalConstraints, color_2, 'LineStyle', 'none');
end

for ii = 1:sls_mpc.T
    color_1 = [178 102 255]./255;
    PolytopePlot.show_convex(tube_sol.tubes{ii}, color_1, 'FaceAlpha', 0.3);
    hold on
end

N = length(tubeTrajSet);
for ii = 1:N
    plot(tubeTrajSet{ii}.x_seq(1,:), tubeTrajSet{ii}.x_seq(2,:), 's-', 'LineWidth', 1.5);
    hold on
end
x0 = sls_mpc.x0;
scatter(x0(1), x0(2),'oy','filled');

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title(['Tube MPC trajectories, T = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

figure;
for ii = 1:N
    plot(tubeTrajSet{ii}.u_seq, 's-', 'LineWidth', 1.5);
    hold on
end
xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
title(['Tube MPC control inputs, T = ',num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);


end

