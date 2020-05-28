function [solCell, sets] = VerifyFeasibility_Coarse(MPC, x0Set, options)
%VERIFYFEASIBILITY_COARSE % grid search over the initial conditions x0 in 
% the set x0Set and check the feasibility of robust coarse SLS MPC on each 
% x0. Each x0 has three status: feasible, infeasible and unverified.

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
    
    [bounds,isFeasible] = MPC.CoarseBisectParams(range, bisect_opt);

    sol.bounds = bounds;
    fprintf('Searching for feasible parameters...\n');
    if isFeasible
        sol.feasibleBounds = 1;
        [feasibleParams] = MPC.CoarseGridSearchParams(bounds, gridDim, gridparams_opt)
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

sets.feasibleSet = feasibleSet; 
sets.infeasibleSet = infeasibleSet;
sets.unverifiedSet = unverifiedSet;


end

