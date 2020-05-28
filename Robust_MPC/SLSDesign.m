classdef SLSDesign < handle
    %SLSDESIGN class of a robust MPC problem equipped with robust SLS,
    %tube MPC and dynamic programming (DP) methods.
    %   Task: solve a robust optimal control problem (OCP) in each MPC loop
    %   Consider robust closed-loop OCP of the uncertain system with given
    %   x(0) and horizon T:
    %   min_pi J_nominal(x(0), pi)
    %   s.t. for all norm bounded disturbance and model uncertainty
    %        norm(w(k)) <= sigma_w, norm(Delta_A)<= eps_A,
    %        norm(Delta_B)<= eps_B:
    %           x(k+1) = (Ahat + Delta_A)x(k) + (Bhat + Delta_B)u(k) + w(k)
    %           u(k) = pi(x(0), x(1), ..., x(k))
    %           x(k) \in X, u(k) \in U, x(T) \in X_T, forall k <= T
    %   For details of the formulation see the reference paper "Robust
    %   closed-loop model predictive control via system level synthesis" by
    %   S.Chen, H.Wang, M.Morari, V.M.Preciado, N.Matni, https://arxiv.org/pdf/1911.06842.pdf

    properties (SetAccess = public)
        % parameters that describe the robust MPC problem
        Ahat; Bhat; % nominal model
        AhatSeq = []; BhatSeq = []; % initialize if the nominal model is time variant; otherwise left empty
		Q; R;  % Q: state weight matrix, R: input weight matrix
        terminalCost; % terminal weight matrix
        T; % horizon
		stateDim; inputDim; % state and input dimensions
        
        % the following state, input and terminal constraints are
        % represented by the Polyhedron class available in MPT3 toolbox
		stateConstraints; inputConstraints;
        disturbanceConstraints;
        terminalConstraints; 
        
		x0; % initial condition
		sigmaW; % norm bound on disturbance w(k)
		noiseNormType = 2; 
        FIR = true; 
		
		options; % a struct that is used for passing parameters
		epsA; epsB; eps; % norm bound on model uncertainty
        
        Phi_x = []; Phi_u =[]; % save the Phi_x, Phi_u solutions after running SLS
        sequence = []; % save state and input sequences
        hyperParamBounds = [];
    end
    
    %% methods section
    methods (Access = public)
        
        %% initialization
		function obj = SLSDesign(params, options)
            % Extracting parameters
            if nargin < 2
                options = struct;
            end
            % nominal system
            Ahat = params.Ahat; Bhat = params.Bhat; 
            x0 = params.x0;
            Q = params.Q; R = params.R; terminalCost = params.terminalCost; 
            T = params.horizon; 
            % level of uncertainty
            eps = params.eps; sigmaW = params.sigmaW;
            % state, input, terminal and disturbance constraints
            stateConstraints = params.stateConstraints; 
            inputConstraints = params.inputConstraints;
            disturbanceConstraints = params.disturbanceConstraints; 
            terminalConstraints = params.terminalConstraints;
            % the norm type of the disturbance
            if ~isfield(params, 'noiseNormType')
                noiseNormType = 2;
            else
                noiseNormType = params.noiseNormType;
            end
            
            % Assign the values of properties
			obj.Ahat = Ahat; 
			obj.Bhat = Bhat;  
            obj.T = T;
            if isfield(params, 'AhatSeq')
                obj.AhatSeq = params.AhatSeq;
                if length(obj.AhatSeq) ~= T
                    warning('Length of sequence not match!\n');
                end
            end
            if isfield(params, 'BhatSeq')
                obj.BhatSeq = params.BhatSeq;
                if length(obj.AhatSeq) ~= T
                    warning('Length of sequence not match!\n');
                end
            end
            
			obj.Q = Q; obj.R = R; 
            if length(eps) == 1
                obj.eps = eps; obj.epsA = []; obj.epsB = [];
            elseif length(eps) == 2
                obj.epsA = eps(1); obj.epsB = eps(2); obj.eps = eps(1) + eps(2);
            end
  			obj.x0 = x0;
			obj.sigmaW = sigmaW; 
            obj.noiseNormType = noiseNormType;
			obj.stateConstraints = stateConstraints;
			obj.inputConstraints = inputConstraints;
            obj.disturbanceConstraints = disturbanceConstraints;
			obj.terminalConstraints = terminalConstraints;
			obj.terminalCost = terminalCost;
			obj.options = options;
			obj.stateDim = size(Ahat, 1);
			obj.inputDim = size(Bhat, 2);
        end
        
        %% solve robust MPC through SLS
        function [sol] = SolveSLSMPC(obj, params, opt) 
            % Inputs:
            %   params: struct of hyperparameters. Fields: .beta, .gamma,
            %   .tau, .alpha (optional)
            %   opt: struct, fields: .verbose, .keepGamma, .keepAffine,
            %   .keepTau, .keepRobustPerf, .keepSafety, .keepObj
            
             if nargin < 3
                opt.verbose = 2;
%                 opt.blockSafetyConstr = 0;
%                 opt.blockBoxConstr = 0;
                opt.keepGamma = 1;
                opt.keepAffine = 1;
                opt.keepTau = 1;
                opt.keepRobustPerf = 1;
                opt.keepSafety = 1;
                opt.keepObj = 1;
             end
               
            % extraction hyperparameters
            tau = params.tau; beta = params.beta; gamma = params.gamma;
            if ~isfield(params, 'alpha')
                alpha = 0.5;
            else
                alpha = params.alpha; 
            end
            
            if ~isfield(opt,'verbose')
                 opt.verbose = 2;
            end 
            verbose = opt.verbose;

        	Ahat = obj.Ahat; Bhat = obj.Bhat;
        	n = obj.stateDim; p = obj.inputDim; T = obj.T; 
			x0 = obj.x0; sigmaW = obj.sigmaW; 
            epsA = obj.epsA; epsB = obj.epsB;
            eps = obj.eps;
       
            % construct SLS 
			sqrtQ = sqrtm(obj.Q);
			sqrtR = sqrtm(obj.R);

			Phi_x = sdpvar( (T + 1) * n, (T + 1) * n, 'full');
			Phi_u = sdpvar( (T + 1) * p, (T + 1) * n, 'full');
            % extract matrix blocks
            Phi = [Phi_x; Phi_u];
            Phi_x_0 = Phi_x(:,1:n);
            Phi_u_0 = Phi_u(:,1:n);
            Phi_x_w = Phi_x(:,n+1:end);
            Phi_u_w = Phi_u(:,n+1:end);
            Phi_0 = Phi(:,1:n);
            Phi_w = Phi(:,n+1:end);
            
            % set of constraints
			constr = [];

			% Construct the objective function using nominal cost
			if ~isempty(obj.terminalCost)
				P = obj.terminalCost;
				sqrtP = sqrtm(P);
			else
				sqrtP = sqrtQ;
			end

			Qset = cell(T + 1, 1); Rset = cell(T + 1, 1);
			Qset{1} = zeros(size(sqrtQ));
			for i = 2 : T 
				Qset{i} = sqrtQ;
			end
			Qset{T + 1} = sqrtP;

			for i = 1:T
				Rset{i} = sqrtR;
            end
            % no penalty on u_T
			Rset{T + 1} = zeros(size(sqrtR));
            
			cost = 0;
            if opt.keepObj
			for i  = 1 : T+1
				cost = cost + norm(Qset{i}*Phi_x((i - 1)*n + 1: i*n, 1:n )*x0, 2)^2 ...
				       + norm(Rset{i}*Phi_u((i - 1)*p + 1: i*p, 1:n )*x0, 2)^2;      
            end
          
            end

			% Structure constraint of Phi_x and Phi_u
			constr = [constr, Phi_x(1:n, 1:n) == eye(n)];

			for k = 1 : T
				constr = [constr, Phi_x( (k - 1)*n + 1: k*n, k*n + 1: end) == zeros(n, (T + 1 - k)*n)];
			end

			for k = 1 : T
				constr = [constr, Phi_u( (k - 1)*p + 1: k*p, k*n + 1: end) == zeros(p, (T + 1 - k)*n)];
            end
            
			% Affine subspace constraint
            if ~isempty(obj.AhatSeq)
               AhatSeq = obj.AhatSeq;
               blkA = blkdiag(AhatSeq{:}, zeros(size(Ahat)));
            else
               blkA = blkdiag(kron(eye(T), Ahat), zeros(size(Ahat)));
            end
            
            if ~isempty(obj.BhatSeq)
               BhatSeq = obj.BhatSeq;
               blkB = blkdiag(BhatSeq{:}, zeros(size(Bhat)));
            else
               blkB = blkdiag(kron(eye(T), Bhat), zeros(size(Bhat)));
            end

			Z = kron(diag(ones(1,T),-1), eye(n));
			
			Id = eye((T + 1)*n);
            
            if opt.keepAffine
                constr = [constr, (Id - Z*blkA)*Phi_x - Z*blkB*Phi_u == Id];
            end
            
			% State and input constraints

            % % construct F
			if ~isempty(obj.stateConstraints)
				Fx = obj.stateConstraints.A; bx = obj.stateConstraints.b;
				nFx = size(Fx, 1); nbx = length(bx);
            else
                warning('must have state constraints');
            end
            
            if ~isempty(obj.terminalConstraints)
				Ft = obj.terminalConstraints.A; bt = obj.terminalConstraints.b;
				nFt = size(Ft, 1); nbt = length(bt);
            else 
                Ft = Fx; bt = bx; nFt = nFx; nbt = nbx;
            end
            
            if ~isempty(obj.inputConstraints)
				Fu = obj.inputConstraints.A; bu = obj.inputConstraints.b;
				nFu = size(Fu, 1); nbu = length(bu);
            else
                warning('must have input constraints');
            end
            
            F = blkdiag(Fx);
            for j = 1:T-1
                F = blkdiag(F, Fx);
            end
            F = blkdiag(F, Ft);
            for j = 1:T
                F = blkdiag(F, Fu);
            end
            % concatenate all input, state and terminal constraints.
            F = blkdiag(F, zeros(nFu, p));
            b = [kron(ones(T,1), bx); kron(ones(1,1), bt); kron(ones(T,1),bu); kron(ones(1,1), zeros(size(bu)))];
            
            nFconstr = size(F,1);
            assert(nFconstr == T*nFx + nFt + (T+1)*nFu);
            assert(size(b,1) == T*nFx + nFt + (T+1)*nFu);
            assert(size(b,2) == 1);          
           
            % construct slack variable
            for j = 1:nFconstr - nFu % no constraint on u_T
                if opt.keepSafety
                   constr = [constr, F(j,:)*Phi_0*x0 + norm(F(j,:)*Phi_w,1)*(1 - tau^T)/(1 - tau)*gamma + beta*sigmaW  <= b(j)];
                end

            if sigmaW ~= 0
                if opt.keepRobustPerf
                    constr = [constr, norm(F(j,:)*Phi_w, inf) + beta*norm(eps*Phi_w, inf) <= beta];
                end
            end
            end

            if isempty(epsA)
                if opt.keepTau
                    constr = [constr, norm(Phi_w, inf)*eps <= tau];
                end
                if opt.keepGamma  
                    constr = [constr, eps*norm(Phi_0*x0, inf) <= gamma];
                end
            else
                if opt.keepTau
                    Phi_temp_1 = [epsA/alpha*Phi_x_w; epsB/(1 - alpha)*Phi_u_w];
                    constr = [constr, norm(Phi_temp_1, inf) <= tau];
                end
                if opt.keepGamma
                    Phi_temp_2 = [epsA/alpha*Phi_x_0; epsB/(1 - alpha)*Phi_u_0];
                    constr = [constr, norm(Phi_temp_2*x0, inf) <= gamma];
                end
            end 

            % solve the problem
      
			ops = sdpsettings('verbose',verbose, 'solver', 'mosek');
            solution = optimize(constr, cost, ops);      
      
			if solution.problem ~= 0
				warning('Numerical error detected. \n');
			end

			sol = struct;
			Phi_x_val = value(Phi_x); Phi_u_val = value(Phi_u);
            fval = value(cost);

            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            sol.fval = fval; 
            sol.tau = tau; sol.alpha = alpha; sol.gamma = gamma; sol.beta = beta;
            sol.status = solution.problem;
            sol.time = solution.yalmiptime;
            sol.solution = solution;

            obj.Phi_x = Phi_x_val; obj.Phi_u = Phi_u_val;  
            yalmip('clear'); % clear all yalmip variables
            
        end
        
        %% bisect to find lower and upper bounds of the hyperparameters gamma, beta, tau
        function [bounds, isFeasible, solvertime] = BisectParams(obj, range, options)
            % Inputs:
            %   range: struct, for the form range.*.lb and range.*.ub where
            %          * = gamma, beta, tau
            %   options: struct, fields: options.init, options.tol,
            %   options.verbose.
            if ~isfield(options, 'init')
                init = struct;
            else
                init = options.init;
            end
            
            if ~isfield(options, 'tol')
                tol = 1e-3;
            else 
                tol = options.tol;
            end
            
            if ~isfield(options, 'verbose')
                verbose = 0;
            else
                verbose = options.verbose;
            end
            
            % hyperparams for obj.SolveSLSMPC
            params = struct;
            params.alpha = 0.5;
            params.beta = 1;
            params.gamma = 0.2;
            params.tau = 0.2;

            isFeasible = 1;            
            bounds = struct;
            solvertime = 0;
            
            % find lower bounds for beta, gamma, tau first
            % Find lower bound on gamma
            fprintf('Finding lower bound on gamma...\n');

            if isfield(init,'gamma') & isfield(init.gamma, 'lb') & ~isempty(init.gamma.lb)
                bounds.gamma.lb = init.gamma.lb;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 1;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.gamma.lb; ub = range.gamma.ub;
                mid = (lb + ub)/2;

                while ub - lb > tol
                    params.gamma = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.gamma.lb = mid;
            end

            % Find lower bound on beta
            fprintf('Finding lower bound on beta...\n');

            if isfield(init,'beta') & isfield(init.beta, 'lb') & ~isempty(init.beta.lb)
                bounds.beta.lb = init.beta.lb;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 1;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.beta.lb; ub = range.beta.ub;
                mid = (lb + ub)/2;
                while ub - lb > tol
                    params.beta = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.beta.lb = mid;
            end
            
            % Find lower bound on tau
            fprintf('Finding lower bound on tau...\n');

            if isfield(init,'tau') & isfield(init.tau, 'lb') & ~isempty(init.tau.lb)
                bounds.tau.lb = init.tau.lb;  
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 1;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;

                while ub - lb > tol
                    params.tau = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.tau.lb = mid;      
            end
            
            % find upperbounds for beta, gamma, tau
            % Find upper bound on beta
            fprintf('Finding upper bound on beta...\n');

            if isfield(init,'beta') & isfield(init.beta, 'ub') & ~isempty(init.beta.ub)
                bounds.beta.ub = init.beta.ub;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 1;
                opt.keepObj = 0;
            
                lb = range.beta.lb; ub = range.beta.ub;
                mid = (lb + ub)/2;
                params.tau = bounds.tau.lb; params.gamma = bounds.gamma.lb;
                while ub - lb > tol
                    params.beta = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        ub = mid; 
                        mid = (lb + ub)/2;
                    else
                        lb = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.beta.ub = mid; 
            end
            
            if bounds.beta.ub < bounds.beta.lb
                isFeasible = 0;
                return
            end
                       
            % find upperbound on tau
            fprintf('Finding upper bound on tau...\n');

            if isfield(init,'tau') & isfield(init.tau, 'ub') & ~isempty(init.tau.ub)
                bounds.tau.ub = init.tau.ub;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 1;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;
                params.gamma = bounds.gamma.lb; params.beta = bounds.beta.lb;
                while ub - lb > tol
                    params.tau = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        ub = mid; 
                        mid = (lb + ub)/2;
                    else
                        lb = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.tau.ub = mid;   
            end
            
            if bounds.tau.ub < bounds.tau.lb
                isFeasible = 0;
                return
            end
            
            % find upperbound on gamma
            fprintf('Finding upper bound on gamma...\n');

            if isfield(init,'gamma') & isfield(init.gamma, 'ub') & ~isempty(init.gamma.ub)
                bounds.gamma.ub = init.gamma.ub;
            else
                tempup = (1 - bounds.tau.ub^obj.T)/(1 - bounds.tau.ub);
                templow = (1 - bounds.tau.lb^obj.T)/(1 - bounds.tau.lb);
                bounds.gamma.ub = tempup*bounds.gamma.lb/templow;
            end
            
            if bounds.gamma.ub < bounds.gamma.lb
                isFeasible = 0;
                return
            end
            fprintf('Bisection finished.\n');
        end
        
        %% grid search for a feasible set of hyperparameters
        function [feasibleParams, solvertime] = GridSearchParams(obj, bounds, gridDim, options)
            % Grid search for a feasible set of hyperparameters (gamma,
            % beta, tau).
            % Inputs:
            %   bounds: struct, fields: bounds.*.lb, bounds.*.ub, with * = beta, tau,
            %   gamma.
            %   gridDim: struct, fields: gridDim.num_* with * = beta, tau, gamma.
            %   options: struct, fields: options.verbose
            
            if ~isfield(options, 'verbose')
                verbose = 0;
            else
                verbose = options.verbose;
            end
                   
            feasibleParams = [];
            solvertime = 0;
            
            solveOpt = struct;
            solveOpt.verbose = verbose;
            solveOpt.keepGamma = 1;
            solveOpt.keepAffine = 1;
            solveOpt.keepTau = 1;
            solveOpt.keepRobustPerf = 1;
            solveOpt.keepSafety = 1;
            solveOpt.keepObj = 0;

            params.alpha = 0.5;
            
            num_tau = gridDim.num_tau; 
            num_gamma = gridDim.num_gamma; 
            num_beta = gridDim.num_beta;
            
            betaRange = linspace(bounds.beta.lb, bounds.beta.ub, num_beta + 2);
            gammaRange = linspace(bounds.gamma.lb, bounds.gamma.ub, num_gamma + 2);
            tauRange = linspace(bounds.tau.lb, bounds.tau.ub, num_tau + 2);
            
            % discard samples right on the boundary (they tend to be
            % infeasible)
            for beta = betaRange(2:end-1)
                for gamma = gammaRange(2:end-1)
                    for tau = tauRange(2:end-1)
                        params.beta = beta; params.gamma = gamma; params.tau = tau;
         
                        fprintf('grid search for SLS MPC feasible params: %f / %d, %f / %d, %f / %d \n', ...
                        find(betaRange == beta), length(betaRange)-2, find(gammaRange == gamma), length(gammaRange)-2, ...
                        find(tauRange == tau), length(tauRange)-2);
                        [sol] = obj.SolveSLSMPC(params, solveOpt); 
                        solvertime = solvertime + sol.solution.solvertime;
                        % return on finding the first feasible set
                        if sol.status == 0
                            feasibleParams = params;
                            fprintf('Feasible parameters found!\n');
                            return
                        end
                    end
                end
            end 
            fprintf('Grid search terminated. Good luck next time! \n');
        end
        
        %% generate trajectories
        function [trajSet, nominalTraj] = SimulateTraj(obj, num_traj, options)
           % simulate closed-loop trajectories under the synthesized controller.
           % num_traj: number of experiments.
           if nargin < 3
               options = struct;
           end
           
           % ismemory = 1: apply model uncertainty with memory
           % otherwise apply memoryless model uncertainty
           if isfield(options, 'ismemory')
               ismemory = options.ismemory;
           else
               ismemory = 1;
           end
           
           T = obj.T; x0 = obj.x0;
           n = obj.stateDim; p = obj.inputDim;
           sigmaW = obj.sigmaW;
           
           Ahat = obj.Ahat; Bhat = obj.Bhat;
           if ~isempty(obj.AhatSeq)
               AhatSeq = obj.AhatSeq;
           else
               AhatSeq = mat2cell(repmat(Ahat, 1, T), size(Ahat,1), ones(1,T)*size(Ahat, 2));
           end
           
           if ~isempty(obj.BhatSeq)
               BhatSeq = obj.BhatSeq;
           else
               BhatSeq = mat2cell(repmat(Bhat, 1, T), size(Bhat,1), ones(1,T)*size(Bhat, 2));
           end
           
           if isempty(obj.Phi_x) | isempty(obj.Phi_u)
                error('The MPC problem must be solved first. \n')
           end
           
           if isempty(obj.epsA)
               epsA = obj.eps; epsB = obj.eps;
           else
               epsA = obj.epsA; epsB = obj.epsB;
           end
           
           % if specified, use the designated parameters
           if isfield(options, 'epsA')
               epsA = options.epsA;
           end
           
           if isfield(options, 'epsB')
               epsB = options.epsB;
           end
           Phi_x = obj.Phi_x; Phi_u = obj.Phi_u;
           % synthesize the feedback controller
           K = Phi_u*inv(Phi_x);
           
           trajSet = cell(1, num_traj);
           for ii = 1:num_traj
               traj = struct;
               x_seq = [x0]; u_seq = []; w_seq = [x0];
               DeltaA_seq = cell(1, T); DeltaB_seq = cell(1, T);
               xsls = x0;
               for jj = 1:T
                   w = rand(n,1).*2*sigmaW - sigmaW;           
                   w_seq = [w_seq w];

                   usls = zeros(p, 1);
                   for kk = 1:jj
                       usls = usls + K((jj - 1)*p + 1:jj*p,(kk - 1)*n + 1:kk*n)*x_seq(:,kk);
                   end
                   u_seq = [u_seq usls];
                    
                   % obtain the next state
                   [Delta_A] = DeltaOperator(n, n, jj, epsA, ismemory);
                   [Delta_B] = DeltaOperator(n, p, jj, epsB, ismemory);
                   DeltaA_seq{jj} = Delta_A;
                   DeltaB_seq{jj} = Delta_B;
                   
                   dist_A = Delta_A*x_seq(:);
                   dist_B = Delta_B*u_seq(:);
                   
                   xsls = AhatSeq{jj}*xsls + dist_A + BhatSeq{jj}*usls + dist_B + w;
                   x_seq = [x_seq xsls];
               end
               traj.x_seq = x_seq;
               traj.u_seq = u_seq;
               traj.w_seq = w_seq;
               traj.DeltaA_seq = DeltaA_seq;
               traj.DeltaB_seq = DeltaB_seq; 
               
               trajSet{ii} = traj;
           end
           
           nominalTraj = struct;
           x_nominal_seq = [x0]; u_nominal_seq = [];
           for ii = 1:T
               u_nominal_seq = [u_nominal_seq Phi_u((ii - 1)*p + 1:ii*p, 1:n)*x0];
               x_nominal_seq = [x_nominal_seq Phi_x((ii)*n + 1:(ii+1)*n, 1:n)*x0];
           end
           nominalTraj.u_nominal_seq = u_nominal_seq;
           nominalTraj.x_nominal_seq = x_nominal_seq;
        end
        
        %% dynamic programming approach
        function [output] = SolveDP(obj, normType)
            % solve the robust OCP with DP. Yalmip toolbox required.
            % see https://yalmip.github.io/example/robustmpc/ for details
            
            if nargin < 2
                warning('Please specify the type of the norm.\n');
                 normType = 1;
            end
            
            Ahat = obj.Ahat; Bhat = obj.Bhat;
            nx = obj.stateDim; nu = obj.inputDim;
            epsA = obj.epsA; epsB = obj.epsB; sigmaW = obj.sigmaW;
            T = obj.T;
            
            Xc = obj.stateConstraints; 
            Uc = obj.inputConstraints;
            Wc = obj.disturbanceConstraints;
            Xf = obj.terminalConstraints;
            
            x = sdpvar(repmat(nx, 1, T), repmat(1, 1, T));
            u = sdpvar(repmat(nu, 1, T), repmat(1, 1, T));
            Delta_A = sdpvar(repmat(nx, 1, T), repmat(nx, 1, T));
            Delta_B = sdpvar(repmat(nx, 1, T), repmat(nu, 1, T));
            w = sdpvar(repmat(nx, 1, T), repmat(1, 1, T));
            
            if isempty(obj.terminalCost)
                J{T} = 0;
            elseif normType == 2
                J{T} = x{T}'*obj.terminalCost*x{T};
            else
                J{T} = norm(x{T}, normType); 
            end

            for ii = T-1:-1:1
                fprintf('Step %d out of %d \n', T - ii, T-1);
                % feasible region
                F = [ismember(x{ii}, Xc), ismember(u{ii}, Uc)];
                if ii == T - 1 & ~isempty(Xf)
                    F = [F, ismember(Ahat*x{ii} + Delta_A{ii}*x{ii} + Bhat*u{ii} + Delta_B{ii}*u{ii} + w{ii}, Xf)];
                end
                F = [F, ismember(Ahat*x{ii} + Delta_A{ii}*x{ii} + Bhat*u{ii} + Delta_B{ii}*u{ii} + w{ii}, Xc)];
            

                F = [F, uncertain(w{ii}), uncertain(Delta_A{ii}), uncertain(Delta_B{ii})];
                F = [F, ismember(w{ii}, Wc), norm(Delta_A{ii}, inf) <= epsA, norm(Delta_B{ii}, inf) <= epsB];

                Jw = replace(J{ii+1},x{ii+1}, Ahat*x{ii} + Bhat*u{ii});
                
                switch normType
                    case 2
                        cost = x{ii}'*obj.Q*x{ii} + u{ii}'*obj.R*u{ii} + Jw;
                    case 1
                        cost = obj.Q(1)*norm(x{ii}, 1) + obj.R(1)*norm(u{ii}, 1) + Jw;
                    case Inf
                        cost = obj.Q(1)*norm(x{ii}, Inf) + obj.R(1)*norm(u{ii}, Inf) + Jw;    
                    otherwise
                        warning('p = 1, 2 or inf must be satisfied!\n');
                        cost = 0;
                end
                 % Solve one-step problem
                [sol{ii},diagnost{ii},Uz{ii},J{ii},Optimizer{ii}] = solvemp(F,cost,[],x{ii},u{ii});
            end
            output = struct;
            output.sol = sol;
            output.diagnost = diagnost;
            output.Uz = Uz;
            output.J = J;
            output.Optimizer = Optimizer;
            
        end
       
        %% tube MPC approach
        function [sol] = SolveTubeMPC(obj, Z, options)
            % solve robust OCP through tube MPC
            % details see "Robust model predictive control using tubes" by
            % W.Langson, I.Chryssochoos, S.V.Rakovic, D.Q.Mayne
            % Z: a Polyhderon instance, the shape of the tube
            fprintf('Solving tube MPC started... \n');
            if nargin < 3
                options = struct;
            end
            
            if isfield(options, 'normType')
                normType = options.normType;
            else
                normType = 2;
            end
               
            if isfield(options, 'verbose')
                verbose = options.verbose;
            else
                verbose = 0;
            end
            
            Ahat = obj.Ahat; Bhat = obj.Bhat;
            nx = obj.stateDim; nu = obj.inputDim;
            Q = obj.Q; R = obj.R; P = obj.terminalCost;
            epsA = obj.epsA; epsB = obj.epsB; sigmaW = obj.sigmaW;
            T = obj.T;
            
            x0 = obj.x0; 
            Xc = obj.stateConstraints; 
            Uc = obj.inputConstraints;
            Wc = obj.disturbanceConstraints;
            Xf = obj.terminalConstraints;
            
            Aw = Wc.A; bw = Wc.b;

            % shape of the tubes
            Z_vertices = Z.V;
            assert(size(Z_vertices,2) == nx);
            
            nJ = size(Z_vertices, 1);
            % set of vertices of Z
            Z_Vset = mat2cell(Z_vertices', [size(Z_vertices, 2)], ones(1, nJ));
            
            % find the H representation of Z
            Az = Z.A; bz = Z.b;
            assert(isempty(Z.Ae));
            
            % for computing the Minkowski difference
            nz = size(Az, 1);
            max_W = zeros(nz, 1);
            if sigmaW ~= 0
                opts = optimoptions(@linprog,'Display', 'off');
                for ii = 1:nz
                   [~, fval] = linprog(-Az(ii,:)', Aw, bw,[],[],[],[],opts );
                   max_W(ii) = -fval;
                end
            end
            
            % create optimization variables
            alpha = sdpvar(ones(1,T), ones(1,T));
            zc = sdpvar(repmat(nx, 1, T), repmat(1, 1, T));
            u0 = sdpvar(nu, 1);
            for ii = 1:T
                U{ii} = sdpvar(repmat(nu, 1, nJ), repmat(1, 1, nJ));
                % X{i}{j}: the j-th vertex of the tube X_i
                X{ii} = sdpvar(repmat(nx, 1, nJ), repmat(1, 1, nJ));
            end

            for ii = 1:T
                for jj = 1:nJ
                   X{ii}{jj} = zc{ii} + alpha{ii}*Z_Vset{jj};
                end
            end
            
            % vertices of the inf-norm ball of Delta_A and Delta_B
            [Delta_A_vertices] = FindMatVertices(size(Ahat), epsA);
            [Delta_B_vertices] = FindMatVertices(size(Bhat), epsB);
            num_Delta_A = length(Delta_A_vertices);
            num_Delta_B = length(Delta_B_vertices);
            
            % construct the optimization problem
            fprintf('Constructing constraints ... \n');
            F = [ismember(u0, Uc)];
            % state and input constraints on tubes
            for ii = 1:T
                F = [F, alpha{ii} >= 0];
                for jj = 1:nJ
                   F = [F, ismember(X{ii}{jj}, Xc), ismember(U{ii}{jj}, Uc)];
                   if ii == T & ~isempty(Xf)
                      F = [F, ismember(X{ii}{jj}, Xf)]; 
                   end
                end
            end

            for ii_Delta_A  = 1:num_Delta_A
                for ii_Delta_B = 1:num_Delta_B
                   for ii = 1:T-1
                      for jj = 1:nJ
                        F = [F, Az*(Ahat + Delta_A_vertices{ii_Delta_A})*X{ii}{jj} + Az*(Bhat + Delta_B_vertices{ii_Delta_B})*U{ii}{jj}...
                             - Az*zc{ii+1} + max_W <= alpha{ii+1}*bz];
                      end          
                   end
                   F = [F, Az*(Ahat + Delta_A_vertices{ii_Delta_A})*x0 + Az*(Bhat + Delta_B_vertices{ii_Delta_B})*u0 ...
                           - Az*zc{1} + max_W <= alpha{1}*bz];
                end
            end

            % construct cost function
            switch normType
                case 2
                    cost_fcn = @(x,u) x'*Q*x + u'*R*u;
                    cost_terminal = @(x) x'*P*x;
                case 1
                    cost_fcn = @(x,u) Q(1)*norm(x, 1) + R(1)*norm(u, 1);
                    cost_terminal = @(x) P(1)*norm(x, 1);
                case Inf
                    cost_fcn = @(x,u) Q(1)*norm(x, Inf) + R(1)*norm(u, Inf);
                    cost_terminal = @(x) P(1)*norm(x, Inf);
                otherwise
                    warning('The norm type has to be 1, 2, or Inf.\n');
            end
            
            fprintf('Constructing objective function... \n');
            cost = cost_fcn(x0, u0);
            for ii = 1:T-1
                for jj = 1:nJ
                    cost = cost + cost_fcn(X{ii}{jj}, U{ii}{jj});
                end
            end
            % add terminal cost
            if ~isempty(P)
                for jj = 1:nJ
                   cost = cost + cost_terminal(X{T}{jj}); 
                end
            end
            
            options = sdpsettings('solver','mosek', 'verbose', verbose);
            sol = struct;
            
            fprintf('Solver started...\n');
            diagnostics = optimize(F, cost);
            
            % extract the tubes
            tubes = cell(1, T);
            for ii = 1:T
                tubes{ii} = value(zc{ii}) + value(alpha{ii})*Z;
            end
            
            % extract the vertices
            vertices_set = cell(1, T);
            for ii = 1:T
               vertices = cell(1, nJ);
               for jj = 1:nJ
                  vertices{jj} = value(X{ii}{jj}); 
               end
               vertices_set{ii} = vertices;
            end
            
            % extract the control inputs
            u_tubes = cell(1,T);
            for ii = 1:T
                vertex_input = cell(1, nJ);
                for jj = 1:nJ
                    vertex_input{jj} = value(U{ii}{jj});
                end
               u_tubes{ii} = vertex_input; 
            end
            
            % save the solution
            sol.diagnostics = diagnostics;
            sol.objective = value(cost);
            sol.u0 = value(u0);
            sol.x0 = x0;
            sol.tubes = tubes;   
            sol.u_tubes = u_tubes;
            sol.vertices_set = vertices_set;
            
            yalmip('clear');
        end
        
        %% find robust invariant set
        function [RIS, isConverge] = RobustInvSet(obj, N_iter)
            % use the "preset intersect current set" algorithm to compute
            % the robust invariance set under given disturbance and model
            % uncertainty. 
            % N_iter: the maximum number of iterations
            % RIS: the set given by the algorithm upon termination
            % isConvergence = 1 if algorithm already convergences; = 0
            % otherwise.
            
           Ahat = obj.Ahat; Bhat = obj.Bhat;
           nx = obj.stateDim; nu = obj.inputDim;
           epsA = obj.epsA; epsB = obj.epsB; sigmaW = obj.sigmaW;
           
           Xc = obj.stateConstraints; 
           Uc = obj.inputConstraints;
           Wc = obj.disturbanceConstraints;
           
           % initialize the iteration   
           Xc = polytope(Xc); Uc = polytope(Uc);
           Wc = polytope(Wc); 
           X_current = Xc;
           
           isConverge = 0;
           for jj = 1:N_iter
               x = sdpvar(nx, 1);
               u = sdpvar(nu, 1);
               Delta_A = sdpvar(nx, nx, 'full');
               Delta_B = sdpvar(nx, nu, 'full');
               w = sdpvar(nx, 1);

               F = [ismember(x, Xc), ismember(u, Uc)];
               F = [F, uncertain(w), uncertain(Delta_A), uncertain(Delta_B)];
               F = [F, ismember(w, Wc), norm(Delta_A, inf) <= epsA, norm(Delta_B, inf) <= epsB];
               F = [F, ismember(Ahat*x + Delta_A*x + Bhat*u + Delta_B*u + w, X_current)];

               cost = 0;

               [sol,diagn,Z,Valuefcn,Optimizer] = solvemp(F,cost ,[],x, u);

               solextract = sol{1};
               % extract the preset X_pre
               X_pre = solextract.Pfinal;
               % take the intersection
               X_new = and(X_pre, X_current);

               if X_new == X_current
                   isConverge = 1;
                   break
               end
               X_current = X_new;

               yalmip('clear');
           end
          
           Vtemp = extreme(X_current);
           RIS = Polyhedron(Vtemp);
        end
        
    end
end

