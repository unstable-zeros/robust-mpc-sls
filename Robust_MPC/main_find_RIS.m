%% find the maximum robust invariant set of the uncertain system.
% it seems that the most computationally cost part is removing the
% redundency of the inequalities or vertices of the lifted polyhedron.

load_data = load('sysdata_DI', 'sysdata');
sysdata = load_data.sysdata;

optUncertain = struct;
uncertainSys = UncertainSystem(sysdata, optUncertain);

Xc = sysdata.stateConstraints;
Uc = sysdata.inputConstraints;

Nstep = 10; % maximum number of iteration
options.robust = 1;
options.minVol = 0.5;
[RIS, diagnostic] = uncertainSys.robustInvariantSet(Xc, Uc, Nstep, options)

save('MaxRIS', RIS);


