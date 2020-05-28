function [] = PlotTraj(sls_mpc, trajSet)
%PLOTTRAJ plot the trajectories given in the set trajSet

T = sls_mpc.T;

figure;
color_1 = [178 102 255]./255;
PolytopePlot.show_convex(sls_mpc.stateConstraints, color_1, 'LineStyle', 'none');
hold on
color_2 = [0 153 76]./255;
if ~isempty(sls_mpc.terminalConstraints)
    PolytopePlot.show_convex(sls_mpc.terminalConstraints, color_2, 'LineStyle', 'none');
end

N = length(trajSet);
for ii = 1:N
    plot(trajSet{ii}.x_seq(1,:), trajSet{ii}.x_seq(2,:), 's-', 'LineWidth', 1.5);
    hold on
end
x0 = sls_mpc.x0;
scatter(x0(1), x0(2),'oy','filled');

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title(['SLS MPC trajectories, T = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

figure;
for ii = 1:N
    plot(trajSet{ii}.u_seq, 's-', 'LineWidth', 1.5);
    hold on
end
xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
title(['SLS MPC control inputs, T = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);


end

