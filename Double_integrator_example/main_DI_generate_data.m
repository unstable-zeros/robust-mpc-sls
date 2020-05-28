%% save the uncertain model
sysdata = struct;
sysdata.Ahat = [1 1; 0 1]; 
sysdata.Bhat = [0.5; 1]; 

nx = size(sysdata.Ahat,2); nu = size(sysdata.Bhat,2);
epsA = 0.02; epsB = 0.05; sigmaW = 0.1;

sysdata.eps = [epsA epsB];
sysdata.epsA = epsA; sysdata.epsB = epsB;
sysdata.Q = eye(nx); sysdata.R = 0.1; 
sysdata.horizon = 10;

sysdata.x0 = [-7;-2]; 
sysdata.sigmaW = sigmaW;

Uc_vertices = [-2; 2];
Uc = Polyhedron(Uc_vertices);

Xc_vertices = [1 1; 1 -1; -1 1; -1 -1]*10;
Xc = Polyhedron(Xc_vertices);

Wc_vertices = [1 1; 1 -1; -1 1; -1 -1]*sigmaW;
Wc = Polyhedron(Wc_vertices);

sysdata.stateConstraints = Xc; sysdata.inputConstraints = Uc;
sysdata.disturbanceConstraints = Wc; 
sysdata.terminalConstraints = [];
sysdata.terminalCost = eye(nx);

% load the maximum robust invariant set as the terminal set
% For a different system, leave terminalConstraints empty first. Then run
% main_find_RIS to find the maximum robust invariant set first and plug it
% in here.
RIStemp = load('RobustRIS_2', 'RIS');
RIS = RIStemp.RIS;
sysdata.terminalConstraints = RIS;

save('sysdata_DI', 'sysdata');
