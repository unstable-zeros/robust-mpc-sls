function [solCell, sets] = VerifyFeasibility(MPC, x0Set, options)
% grid search over the initial conditions x0 in the set x0Set and check the
% feasibility of robust SLS MPC on each x0. Each x0 has three status:
% feasible, infeasible and unverified.

sets = struct;
num_sample = size(x0Set,1);
solCell = cell(1,num_sample);
feasibleSet =[]; infeasibleSet = []; unverifiedSet =[]; 

range = options.range; 
bisect_opt = options.bisect_opt;
gridDim = options.gridDim;
gridparams_opt = options.gridparams_opt;

for ii = 1 : num_sample
    fprintf('Loop %d out of %d \n', ii, num_sample);
    tic
    x0 = x0Set(ii,:)';
    MPC.x0 = x0;
    
    sol = struct;
    sol.x0 = x0;
    sol.bounds = [];
    sol.feasibleParams = [];
    sol.feasibleBounds = 0;
    sol.findFeasibleSol = 0;

    fprintf('Bisecting to find bounds...\n');
    
    [bounds,isFeasible] = MPC.BisectParams(range, bisect_opt);

    sol.bounds = bounds;
    fprintf('Searching for feasible parameters...\n');
    if isFeasible
        sol.feasibleBounds = 1;
        [feasibleParams] = MPC.GridSearchParams(bounds, gridDim, gridparams_opt)
        if ~isempty(feasibleParams)
            sol.findFeasibleSol = 1;
            feasibleSet = [feasibleSet; x0'];
        else
            unverifiedSet = [unverifiedSet; x0'];
        end
        sol.feasibleParams = feasibleParams;
    else
        infeasibleSet = [infeasibleSet; x0'];
    end
    solCell{ii} = sol;

    toc
end

sets.feasibleSet = feasibleSet; % set of feasible x0
sets.infeasibleSet = infeasibleSet; % set of infeasible x0
sets.unverifiedSet = unverifiedSet; % set of unverified x0

end

