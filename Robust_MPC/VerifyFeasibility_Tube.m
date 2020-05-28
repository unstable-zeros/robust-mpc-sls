function [solCell, classSet] = VerifyFeasibility_Tube(MPC, x0Set, Z, options)
% grid search over the initial conditions x0 in the set x0Set and check the
% feasibility of tube MPC on each x0. Each x0 has two status:
% feasible and infeasible.

num_sample = size(x0Set,1);
solCell = cell(1,num_sample);
feasibleSet =[]; infeasibleSet = [];
classSet = struct;
for ii = 1 : num_sample
    fprintf('Loop %d out of %d \n', ii, num_sample);
    tic
    x0 = x0Set(ii,:)';
    MPC.x0 = x0;
    
    [sol] = MPC.SolveTubeMPC(Z, options);
    
    diagnost = sol.diagnostics.problem;
    if diagnost == 0
       feasibleSet = [feasibleSet; x0'];
    else
       infeasibleSet = [infeasibleSet; x0']; 
    end
    
    solCell{ii} = sol;
end

classSet.feasibleSet = feasibleSet; % set of feasible x0
classSet.infeasibleSet = infeasibleSet; % set of infeasible x0

end

