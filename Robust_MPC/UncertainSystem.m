classdef UncertainSystem < handle
    %UNCERTAINSYSTEM auxiliary system class that records the dynamics and
    %uncertainty; mainly used for computing robust invariant sets and
    %finding disturbance invariant tubes. Two functions in this class are
    %important: minInvSet() to find the disturbance invariant set and
    %robustInvariantSet() to find the (open-loop) robust invariant set.
    
    properties (SetAccess = public)
        A; B; % nominal model
        Q; R; % Q: state weight matrix; R: input weight matrix
        DA_vertices; DB_vertices; % vertices of model uncertainty Delta_A, Delta_B
        epsA; epsB; options; % norm ball bounds on Delta_A, Delta_B
        W; % the set that bounds disturbance w(k)
        CLRIS; % closed loop robust invariant set
        OLRIS; % open loop robust invariant set
        K; % LQR controller
    end
    
    methods
        function obj = UncertainSystem(params, options)
            %UNCERTAINSYS Construct an uncertain system based on the robust
            %MPC problem data.
            %   CLRIS and OLRIS save the closed-loop robust invariant set
            %   and open-loop RIS separately. 
            
            obj.A = params.Ahat; obj.B = params.Bhat;
            obj.Q = params.Q; obj.R = params.R;
            D_A_dim = size(params.Ahat); D_B_dim = size(params.Bhat);
            obj.DA_vertices = FindMatVertices(D_A_dim, params.epsA);
            obj.DB_vertices = FindMatVertices(D_B_dim, params.epsB);
            obj.epsA = params.epsA; obj.epsB = params.epsB;
            obj.W = params.disturbanceConstraints;
            obj.options = options;
            obj.CLRIS = []; obj.OLRIS = []; 
            obj.K = [];
        end
        
        %% find K through LQR
        function [K, P] = findK(obj, Q, R)
            [K, P] = dlqr(obj.A, obj.B, Q, R);
        end
        
        %% find minimum disturbance invariant set Z_inv
        function [output, isConverge] = minInvSet(obj, N)
            % N: maximum number of iterations
            W0 = obj.W;
            A = obj.A;
            B = obj.B;
            n = size(A, 1); m = size(B, 2);
            [K, ~] = obj.findK(obj.Q, obj.R);
            Acl = A - B*K;
            W = W0;
            isConverge = 0;
            volList = zeros(1, N+1);
            volList(1) = W.volume();
            for ii = 1:N
               W = W + Acl^ii*W0; % find disturbance invariant set using iterations
               W.minHRep(); % Important: call minHRep to reduce redundant inequalities
               volList(ii+1) = W.volume();
               if abs(volList(ii+1) - volList(ii)) <1e-2
                   isConverge = 1;
                   break;
               end
            end
            output = W;
        end
       
        %% preset of autonomous system
        function [PreS] = presetAuto(obj, Xc, sysA)
            F = Xc.A;  f= Xc.b;
            PreS = Polyhedron(F*sysA, f);
        end
        
        %% preset of autonomous system under model uncertainty
        function [PreS] = presetAutoRobust(obj, Xc, options)
            % Enumeration of vertices of the uncertainty set is applied.
            robustXc = Xc - obj.W;
            F = robustXc.A; f = robustXc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;
            
            K = options.K;
            
            n = size(A,2); m = size(B, 2);
            numDA_vertices = length(DA_vertices); numDB_vertices = length(DB_vertices);
            
            Fbar = []; fbar = [];
            for ii = 1:numDA_vertices
                for jj = 1:numDB_vertices
                    constr = [F*((A + DA_vertices{ii})+(B + DB_vertices{jj})*K)];
                    Fbar = [Fbar; constr]; fbar = [fbar; f];
                end
            end
            
            PreS = Polyhedron(Fbar, fbar);
 
        end
        
        %% preset
        function [PreS] = preset(obj, Xc, Uc, options)
           % compute the preset of Xc for all possible u \in Uc, using
           % projection methods
            if nargin < 4
                options = obj.options;
            end
          
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;

            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 

            Fbar = [F*A F*B; zeros(nU, n) G]; fbar = [f; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);

            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims);
        end

        %% presetRobust
        function [PreS] = presetRobust(obj, Xc, Uc, options)
            % compute the preset of Xc for all possible u \in Uc
            if nargin < 4
                options = obj.options;
            end

            robustXc = Xc - obj.W;
            F = robustXc.A; f = robustXc.b;
            
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;

            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 
               
            numDA_vertices = length(DA_vertices); numDB_vertices = length(DB_vertices);
            
            Fbar = []; fbar = [];
             
            for ii = 1:numDA_vertices
                for jj = 1:numDB_vertices
                    constr = [F*(A + DA_vertices{ii}) F*(B + DB_vertices{jj})];
                    Fbar = [Fbar; constr]; fbar = [fbar; f];
                end
            end
                
            newRow = [zeros(nU, n) G];
            Fbar = [Fbar; newRow]; fbar = [fbar; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);
            tic
            fprintf('remove lift poly redundancy.\n');
            liftedPolyhedron.minHRep();      
            toc
            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims); 
            tic
            fprintf('remove preS redundancy.\n');
            PreS.minHRep();
            toc
        end
        
        
        %% robustInvariantset
        function [RIS, diagnostic] = robustInvariantSet(obj, Xinit, Uc, Nstep, options)
            % Xc is the given initial set on states; Uc is the constraint
            % on u; Nstep is the maximum step simulated forward.
            % This function computes the robust control invariant set.
            if nargin < 4
                Nstep = 10;
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            elseif nargin < 5
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            end
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            volList = zeros(1, Nstep+1);
            volList(1) = Xinit.volume();
            for ii = 1:Nstep
                fprintf('Invariant set iter %d/%d \n', ii, Nstep);
                diagnostic.runningStep = ii;
               
                figure; Xinit.plot(); 
                title(['step = ', num2str(ii), ' total step = ', num2str(Nstep)]);
%                 pause(0.5);
                if options.robust == 0
                    preS = obj.preset(Xinit, Uc, options);
                else 
                    preS = obj.presetRobust(Xinit, Uc, options);
                end
                X_new = and(Xinit, preS);
                X_new.minHRep();
             
                if X_new == Xinit
                    diagnostic.converge = 1;
                    break;
                end
                Xinit = X_new;
                volList(ii+1) = Xinit.volume();
                
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    break;
                end
            end
            RIS = Xinit;
            obj.OLRIS = RIS;
        end
        
        %% Robust invariant set for closed-loop systems
        function [RIS, diagnostic] = robustInvariantSetClosedLoop(obj, Xc, Uc, Nstep, options)
            % compute the robust invariant set of the closed-loop system
            % A+BK
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;
            
            A = obj.A; B = obj.B; K = options.K;
            
%             [F_unify, G_unify, ~] = convert_Poly2Mat(Xc, Uc);
%             Xc_new = Polyhedron(F_unify + G_unify*K, ones(size(F_unify, 1), 1));
            Fbar = [F; G*K]; fbar = [f; g];
            Xc_new = Polyhedron(Fbar, fbar);
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            Xinit = Xc_new;
            sysA = A + B*K;
            
            for ii = 1:Nstep
                diagnostic.runningStep = ii;
                if options.plot ~= 0
                    clf;
                    Graphics.show_convex(Xinit, 'r');
                    if isfield(options, 'xrange') && ~isempty(options.xrange)
                         xlim(options.xrange);
                    end
                    if isfield(options, 'yrange') && ~isempty(options.yrange)
                        ylim(options.yrange);
                    end
%                     Xinit.plot(); 
                    title(['step = ', num2str(ii), ' total step = ', num2str(Nstep)]);
                    pause(0.5);
                end
                if options.robust ~= 1
                    PreS = obj.presetAuto(Xinit, sysA);
                else 
                    PreS = obj.presetAutoRobust(Xinit, options);
                end
                X_new = and(Xinit, PreS);
                if X_new == Xinit
                    diagnostic.converge = 1;
                    break;
                end
                Xinit = X_new;
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    break;
                end
            end
            RIS = Xinit;
            obj.CLRIS = RIS;
            obj.K = K;
        end
        
        %% validate the robust invariant set
        function [x_seq] = validateClosedLoopRIS(obj, Nstep, Ntraj, options)
                RIS = obj.CLRIS;
                Graphics.show_convex(RIS, 'r');
                
                sigmaW = options.sigmaW;
                
                K = obj.K;
                if isempty(RIS) | isempty(K)
                    warning('compute RIS first.\n');
                    return
                end
                sampleList = RIS.grid(20);
                Nsample = size(sampleList,1);
                
                DA_vertices = obj.DA_vertices; 
                DB_vertices = obj.DB_vertices;
                hatA = obj.A; hatB = obj.B;
                
                n = size(hatA, 1);
                
                for kk = 1:Ntraj
                index = randi([1 Nsample]);
                xinit = sampleList(index, :)';

                alpha = rand(1, length(DA_vertices));
                alpha = alpha/sum(alpha);
                beta = rand(1, length(DB_vertices));
                beta = beta/sum(beta);
                
                deltaA = zeros(size(hatA)); deltaB = zeros(size(hatB));
                for ii = 1:length(DA_vertices)
                    deltaA = deltaA + alpha(ii)*DA_vertices{ii};
                end
                    
                for ii = 1:length(DB_vertices)
                    deltaB = deltaB + beta(ii)*DB_vertices{ii};
                end
                
                x_seq = [xinit];
                x_now = xinit;
                for jj = 1:Nstep
                    w_now = rand(n,1).*sigmaW - sigmaW;
                    x_next = ((hatA + deltaA) + (hatB + deltaB)*K)*x_now + w_now;
                    x_seq = [x_seq x_next];
                end
                
                Graphics.show_trajectory(x_seq, 'ys-');   
%                 Graphics.show_trajcetory([xinit], 'gs-');
                end
        end
        
        %% find the Chebyshev center and radius of a polytope
        function [x_c, R, sol] = ChebyshevCenter(obj, Xc)
            F = Xc.A; f = Xc.b;
            
            n = size(F,2);  
            r = sdpvar(1,1); 
            x = sdpvar(n, 1, 'full');
            
            numConstr = size(F, 1);
            constr = [];
            for ii = 1: numConstr 
                constr = [constr, F(ii,:)*x + r*norm(F(ii,:)',2) <= f(ii)];
            end
            constr = [constr, r >= 0];
            
            yalmip_options = sdpsettings('solver','mosek');
            obj = -r;
            sol = optimize(constr,obj,yalmip_options);
            x_c = value(x); R = value(r);
        end
        
        %% Find the maximal RIS through MPC
        function [RIS, PWAcontroller] = RISviaMPC(obj, params, opt)
            if ~isfield(opt, 'plot')
                opt.plot = 1;
            end
            if ~isfield(opt, 'Niter')
                Niter = 20;
            else 
                Niter = opt.Niter;
            end 
            if ~isfield(opt,'tuningN')
                opt.tuningN = 2;
            end
            hatA = obj.A; hatB = obj.B;
            nx = size(obj.A,2); nu = size(obj.B, 2);
            Q = params.Q; R = params.R; 
            if ~isfield(params, 'terminalCost') | isempty(params.terminalCost)
               Qf = Q;
            end
            
            N = opt.tuningN;
            
            if N ~= 1
                u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
            else
                utemp = sdpvar(repmat(nu,1,N),repmat(1,1,N));
                u = {utemp};
            end
            
            x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

            Uc = params.inputConstraints;
            Xc = params.stateConstraints;
            Wc = params.disturbanceConstraints;
            Tc = Xc;
            TcRobust = Tc - Wc;
            XcRobust = Xc - Wc;

            epsA = obj.epsA; epsB = obj.epsB;
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;
            numDA = length(DA_vertices); numDB = length(DB_vertices);
            
            L = numDA*numDB;

            if L ~= 1
                xnextRobust = sdpvar(repmat(nx, 1, L), repmat(1, 1, L));
            else
                xnextRobusttemp = sdpvar(repmat(nx, 1, L), repmat(1, 1, L));
                xnextRobust = {xnextRobusttemp};
            end
            
            % convert Polyhedron instances to polytopes
            Xc = polytope(Xc); Uc = polytope(Uc); Wc = polytope(Wc); 
            Tc = polytope(Tc); TcRobust = polytope(TcRobust); XcRobust = polytope(XcRobust);
            
            for tt = 1:Niter
            % construct constraints
            constraints = [];
            objective = 0;
            for k = 1:N
                 objective = objective +x{k}'*Q*x{k} + u{k}'*R*u{k};
                 constraints = [constraints, x{k+1} == hatA*x{k} + hatB*u{k}];
                 constraints = [constraints, ismember(x{k+1}, XcRobust)];
                 constraints = [constraints, ismember(u{k},Uc)];
            end
            constraints = [constraints, ismember(x{1}, Xc)];

            objective = objective + x{N+1}'*Qf*x{N+1};

            % construct robust contraints on x1
            for ii = 1:numDA
                for jj = 1:numDB
                    constraints = [constraints, xnextRobust{(ii - 1)*numDB + jj} == (hatA + DA_vertices{ii})*x{1} + (hatB + DB_vertices{jj})*u{1}];
                end
            end

            for kk = 1:L
                constraints = [constraints, ismember(xnextRobust{kk}, XcRobust), ismember(xnextRobust{kk},TcRobust)];
            end
            ops = sdpsettings('solver','mosek','verbose', 1);
            [sol,diagn,Z,Valuefcn,Optimizer] = solvemp(constraints,objective ,[],x{1},u{1});

            solextract = sol{1};
            preset = solextract.Pfinal;
            TcNew = and(preset, Tc);
            
            if opt.plot == 1
                figure; plot(TcNew);
            end
            
            TcNew == Tc
            if TcNew == Tc
                break;
            end
            % iterations start here
            Tc = TcNew; TcRobust = Tc - Wc;
            end
            
            % extract the solutions
            Vtemp = extreme(Tc);
            RIS = Polyhedron(Vtemp);
            PWAcontroller = Optimizer;
            
            figure; RIS.plot();
            xlabel('$x_1$','Interpreter', 'LaTex', 'Fontsize',18); 
            ylabel('$x_2$','Interpreter', 'LaTex','Fontsize',18);
            title('RCIS with PWA Local Controller','Interpreter', 'LaTex','Fontsize',18);

            figure;plot(Valuefcn);
            xlabel('$x_1$','Interpreter', 'LaTex', 'Fontsize',18); 
            ylabel('$x_2$','Interpreter', 'LaTex','Fontsize',18);
            zlabel('$J^\star(x)$','Interpreter', 'LaTex','Fontsize',18);
            title('Parametric Optimal Objective Functions','Interpreter', 'LaTex','Fontsize',18);

            figure;plot(Optimizer);
            xlabel('$x_1$','Interpreter', 'LaTex', 'Fontsize',24); 
            ylabel('$x_2$','Interpreter', 'LaTex','Fontsize',24);
            zlabel('$u_0(x)$','Interpreter', 'LaTex','Fontsize',24);
            
            obj.OLRIS = RIS;
        end 
       
        
    end
    
end

