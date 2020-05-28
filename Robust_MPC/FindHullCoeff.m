function [coef, isfeasible] = FindHullCoeff(x, vertices, normType)
%FINDHULLCOEFF find the convex combination coefficients of x w.r.t. the
%vertices. Used for implementing tube MPC control inputs.

if nargin < 3
    normType = 2;
end

num_vertices = size(vertices, 1);
n = size(vertices, 2);
assert(size(x, 1) == n);

lambda = sdpvar(num_vertices, 1);

F = [vertices'*lambda == x, lambda >= 0, ones(1, num_vertices)*lambda == 1];
cost = norm(lambda, normType);

ops = sdpsettings('verbose',0, 'solver', 'mosek');
solution = optimize(F, cost, ops);   
coef = value(lambda);
isfeasible = solution.problem;

yalmip('clear');
end

