%% Model setup 
L = 5000; 
lambda = 500; %10 wavelengths

p = 512;

v0 = 2000; 
freq = v0/lambda;

ppw = p/(L/lambda+2);

model.n = round(L/lambda*ppw*ones(1,3));
model.d = lambda/ppw*ones(1,3);
model.o = [0 0 0];
model.unit = 'm/s';
v = v0*ones(model.n);

nx = model.n(1); ny = model.n(2); nz = model.n(3);

%% Set up source/receiver grid
% options for solving the helmholtz equations
lsopts = LinSolveOpts();
lsopts.tol = 1e-6;
lsopts.maxit = 10000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

% pdefunc options
pdeopts = PDEopts();
pdeopts.helm_dt = model.d;
pdeopts.debug_mode = true;
pdeopts.helm_pml_max = 40;
pdeopts.helm_mat_free = true;
pdeopts.helm_scheme = PDEopts.HELM_OPERTO27;
pdeopts.numcompsrc = 1;

opts = struct;
opts.pdefunopts = pdeopts;
opts.lsopts = lsopts;
opts.subsample_model = false;
opts.disp_progress = true;

[x,y,z] = odn2grid(model.o,model.d,model.n);

% Transmission experiment
%model.xsrc = x(round(linspace(3,nx-3,numsrc)));
%model.ysrc = y(round(linspace(3,ny-3,numsrc)));
model.xsrc = x(3); 
model.ysrc = y(3);
model.zsrc = min(z)+model.d(3);
model.xrec = x(3:1:end-2);
model.yrec = y(3:1:end-2);
model.zrec = min(z)+model.d(3);

% source wavelet time shift
model.t0   = 0; 

% ricker wavelet peak frequency
model.f0 = 10;

% frequencies in Hz
model.freq = 4;

nrx = length(model.xrec); nry = length(model.yrec);  
nfreq = length(model.freq);
nsx = length(model.xsrc); nsy = length(model.ysrc);
nsrc = nsx*nsy;
nrec = nrx*nry;

Q = speye(nsrc);   
sz_data = [nrx,nry,nsx,nsy,nfreq];

%% 
opts.dt = model.d;
opts.pml_max = inf;
opts.solve_opts = lsopts;
opts.mat_free = true;
opts.n_threads = 10;
opts.disp_output = true;
opts.scheme = PDEopts.HELM3D_OPERTO27;
[H,comp_grid] = discrete_helmholtz(v,model,freq,opts);

[xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
sinc_window = 4;
Ps_x = opSincInterp(model.xsrc,xt,sinc_window);
Ps_y = opSincInterp(model.ysrc,yt,sinc_window);
Ps_z = opSincInterp(model.zsrc,zt,sinc_window);
Ps = opKron(Ps_z,Ps_y,Ps_x);

q = Ps*( full(Q(:,1)) );

% You need to measure the memory usage of matlab outside of matlab 
tic,u = H\q; disp(toc);

