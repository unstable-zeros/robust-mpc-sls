classdef CoarseSLSDesign < handle
    %UNTITLED Implementation of coarse SLS to solve robust optimal control
    %problems
    %   This class is mainly for the comparison with the robust SLS
    %   approach proposed in "Robust closed-loop model predictive control 
    %   via system level synthesis" (see https://arxiv.org/pdf/1911.06842.pdf). 
    %   A good reference for coarse SLS is
    %   "System level synthesis" by J.Anderson, J.Doyle, S.Low, N.Matni. (see
    %   https://arxiv.org/pdf/1904.01634.pdf)
    
    properties (SetAccess = public)
        % parameters that describes a robust optiaml control problem
        Ahat; Bhat; AhatSeq = []; BhatSeq = []; % nominal model
		Q; R; % Q: state weight matrix, R: input weight matrix
        terminalCost; % terminal weight matrix
        
        T; % horizon
		stateDim; inputDim; % state and input dimensions
        
        % the following state, input and terminal constraints are
        % represented by the Polyhedron class available in MPT3 toolbox
		stateConstraints; inputConstraints;
        disturbanceConstraints;
        terminalConstraints; 
        
		x0;
		sigmaW; % norm bound on disturbance w(k)
		noiseNormType = 2; 
        FIR = true;
		
		options; 
		epsA; epsB; eps; % norm bound on model uncertainty
        Phi_x = []; Phi_u =[]; % save the Phi_x, Phi_u solutions after running coarse SLS
        sequence = [];
        hyperParamBounds = [];
    end
    
    methods (Access = public)
        %% Initialization
        function obj = CoarseSLSDesign(params, options)
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
        
        %% Solve coarse SLS robust optimal control problem
        function [sol] = SolveCoarseSLSMPC(obj, params, opt)
            if nargin < 3 
                opt = struct;
                opt.keepSafety = 1;
                opt.keepTau = 1;
                opt.keepAffine = 1;
                opt.keepObj = 1;
            end
            
            tau = params.tau;
            if isfield(params, 'alpha')
                alpha = params.alpha;
            else
                alpha = 1/2;
            end
            
            if ~isfield(opt, 'keepObj')
                opt.keepObj = 1;
            end
            
            if isfield(opt,'verbose')
               verbose = opt.verbose;
            else
               verbose = 2;
            end
        		
        	Ahat = obj.Ahat; Bhat = obj.Bhat;
        	n = obj.stateDim; p = obj.inputDim; T = obj.T; 
			x0 = obj.x0; sigmaW = obj.sigmaW; 
			noiseNormType = obj.noiseNormType;
            epsA = obj.epsA; epsB = obj.epsB;
            

			sqrtQ = sqrtm(obj.Q);
			sqrtR = sqrtm(obj.R);

			Phi_x = sdpvar( (T + 1) * n, (T + 1) * n, 'full');
			Phi_u = sdpvar( (T + 1) * p, (T + 1) * n, 'full');

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
			Z = kron(diag(ones(1,T),-1), eye(n));
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

			Id = eye((T + 1)*n);
            if opt.keepAffine
                constr = [constr, (Id - Z*blkA)*Phi_x - Z*blkB*Phi_u == Id];
            end
            
			% State and input constraints
			if noiseNormType == 0
				normType = 1;
			else
				normType = noiseNormType/(noiseNormType - 1);
            end
            
            
            if opt.keepSafety
			if ~isempty(obj.stateConstraints)
				Fx = obj.stateConstraints.A; bx = obj.stateConstraints.b;

				nFx = size(Fx, 1); nbx = length(bx);

				% Assume Fx*x0 <= bx
				% State constraints

				for i = 2 : T + 1 % x0 is not considered and assumed to be feasible
					for j = 1 : nFx
						add_1 = Fx(j,:)*Phi_x((i - 1)*n + 1: i*n, 1:n )*x0 + ...
								norm(Fx(j,:)*Phi_x( (i - 1)*n + 1:i*n, n + 1:i*n), 1)*sigmaW;
						add_2 = norm(Fx(j,:)*Phi_x( (i - 1)*n + 1:i*n, 1:i*n ) ,1)*tau*(1 - tau^T)/(1-tau)*max(sigmaW, norm(x0,inf));
						constr = [constr, add_1 + add_2 <= bx(j)];
					end
                end
            end

            if ~isempty(obj.inputConstraints)
				Fu = obj.inputConstraints.A; bu = obj.inputConstraints.b;

				nFu = size(Fu, 1); nbu = length(bu);

				% Assume Fx*x0 <= bx
				% State constraints

				for i = 1 : T  % x0 is not considered and assumed to be feasible
					for j = 1 : nFu
						add_1 = Fu(j,:)*Phi_u((i - 1)*p + 1: i*p, 1:n )*x0 + ...
								norm(Fu(j,:)*Phi_u( (i - 1)*p + 1:i*p, n + 1:i*n), 1)*sigmaW;
						add_2 = norm(Fu(j,:)*Phi_u( (i - 1)*p + 1:i*p, 1:i*n ) ,1)*tau*(1 - tau^T)/(1-tau)*max(sigmaW, norm(x0,inf));
						constr = [constr, add_1 + add_2 <= bu(j)];
					end
                end
            end

			if ~isempty(obj.terminalConstraints)
				Ft = obj.terminalConstraints.A; bt = obj.terminalConstraints.b;

				nFt = size(Ft, 1); nbt = length(bt);

				% Assume Fx*x0 <= bx
				% State constraints

				for j = 1 : nFt
					add_1 = Ft(j,:)*Phi_x(T*n + 1: (T + 1)*n, 1:n )*x0 + ...
							norm(Ft(j,:)*Phi_x( T*n + 1: (T + 1)*n, n + 1:(T + 1)*n), 1)*sigmaW;
					add_2 = norm(Ft(j,:)*Phi_x( T*n + 1:(T + 1)*n, 1:(T + 1)*n ) ,1)*tau*(1 - tau^T)/(1-tau)*max(sigmaW, norm(x0,inf));
					constr = [constr, add_1 + add_2 <= bt(j)];
				end
            end
            end
            
            
            if opt.keepTau
                Phixu = [ epsA/alpha*Phi_x; epsB/(1 - alpha)*Phi_u];
                constr = [constr, norm(Phixu, inf) <= tau];
            end

			ops = sdpsettings('verbose',verbose, 'solver', 'mosek');
			solution = optimize(constr, cost, ops);      

			if solution.problem ~= 0
				warning('Numerical error detected. \n');
			end

			sol = struct;
			Phi_x_val = value(Phi_x); Phi_u_val = value(Phi_u);
            fval = value(cost);

            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            sol.fval = fval; sol.tau = tau; sol.alpha = alpha;
            sol.status = solution.problem;

            obj.Phi_x = Phi_x_val; obj.Phi_u = Phi_u_val;
            
            yalmip('clear');
        end
        
        %% Bisect to find feasible tau
        function [bounds, isFeasible] = CoarseBisectParams(obj, range, options)
            % Inputs:
            %   range: struct, of the form range.tau.lb and range.tau.ub 
            %   options: struct, fields: options.init, options.tol,
            %   options.verbose.
            %   There is only one hyperparameter tau to bisect.
                        
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
            params.tau = 1;
            
            isFeasible = 1;
            bounds = struct;
            
            % find lower bound on tau
            fprintf('Finding lower bound on tau...\n');
            if isfield(init,'tau') & isfield(init.tau, 'lb') & ~isempty(init.tau.lb)
                bounds.tau.lb = init.tau.lb;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepAffine = 1;
                opt.keepTau = 1;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;

                while ub - lb > tol
                    params.tau = mid;

                    sol = obj.SolveCoarseSLSMPC(params, opt);

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
            
            % find upperbound on tau
            fprintf('Finding upperbound on tau...\n');
            if isfield(init,'tau') & isfield(init.tau, 'ub') & ~isempty(init.tau.ub)
                bounds.tau.ub = init.tau.ub;
            else
                
            opt = struct;
            opt.verbose = verbose;
            opt.keepAffine = 1;
            opt.keepTau = 0;
            opt.keepSafety = 1;
            opt.keepObj = 0;
            
            lb = range.tau.lb; ub = range.tau.ub;
            mid = (lb + ub)/2;
            
            while ub - lb > tol
                params.tau = mid;

                sol = obj.SolveCoarseSLSMPC(params, opt);
  
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
            end
            fprintf('Bisection finished.\n');
        end
        
        %% Grid search for feasible tau
        function [feasibleParams] = CoarseGridSearchParams(obj, bounds, gridDim, options)
            % grid search between the lower and upper bounds to find a
            % feasible hyperparameter tau.
            
            if ~isfield(options, 'verbose')
                verbose = 0;
            else
                verbose = options.verbose;
            end

            feasibleParams = [];
            
            num_tau = gridDim.num_tau;
            if bounds.tau.ub < bounds.tau.lb 
                feasibleParams = [];
                return
            end
            
            solveOpt = struct;
            solveOpt.keepAffine = 1;
            solveOpt.keepTau = 1;
            solveOpt.keepSafety = 1;
            solveOpt.keepObj = 0;
            solveOpt.verbose = verbose;
            
            params.alpha = 1/2;
            
            tauRange = linspace(bounds.tau.lb, bounds.tau.ub, num_tau+2);
            for ii = tauRange(2:end-1)
                fprintf('grid search for SLS-vanilla feasible params: %d / %d \n', ii, length(tauRange)-2);
                params.tau = ii;
                sol = obj.SolveCoarseSLSMPC(params, solveOpt);
                if sol.status == 0
                    feasibleParams = params;
                    fprintf('Feasible parameters found!\n');
                    return
                end
            end
            fprintf('Grid search terminated. Good luck next time! \n');
        end
    
    end
end

