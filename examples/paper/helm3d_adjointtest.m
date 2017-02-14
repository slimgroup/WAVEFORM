% Adjoint test for FWI framework

%% Model setup 
p = 50; 
seed = 1843784;
model.n = [p p p];
model.d = [25 25 25];
model.o = [0 0 0];
model.unit = 'm/s';
v = 2000*ones(model.n);

nx = model.n(1); ny = model.n(2); nz = model.n(3);
rng(seed);

%% Set up source/receiver grid
% options for solving the helmholtz equations
lsopts = LinSolveOpts();
lsopts.tol = 1e-10;
lsopts.maxit = 10000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

% pdefunc options
pdeopts = PDEopts();
pdeopts.helm_dt = model.d;
pdeopts.debug_mode = false;
pdeopts.helm_pml_max = 10;
pdeopts.helm_mat_free = true;
pdeopts.helm_scheme = PDEopts.HELM_OPERTO27;
pdeopts.numcompsrc = 1;

opts = struct;
opts.pdefunopts = pdeopts;
opts.lsopts = lsopts;
opts.subsample_model = false;
opts.disp_progress = false;

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
freq = 4;
opts.dt = model.d;
opts.pml_max = 10;
opts.pml = 10;
opts.solve_opts = lsopts;
opts.mat_free = true;
opts.n_threads = 10;
opts.disp_output = false;
opts.scheme = PDEopts.HELM3D_OPERTO27;
[H,comp_grid] = discrete_helmholtz(v,model,freq,opts);

x = randn(size(H,2),1) + 1i*randn(size(H,2),1);
y = randn(size(H,2),1) + 1i*randn(size(H,2),1);

s = H*x; s = y'*s;
t = H'*y; t = x'*t;
t = conj(t);
disp('Helmholtz matrix');
disp('<H*x,y>');
disp(num2str(s,'%3.15e'));
disp('<x,H''*y>');
disp(num2str(t,'%3.15e'));
disp('Relative difference');

rel_diff = abs(s-t)/max(abs(s),abs(t));
disp(rel_diff);

D = F3d(v,Q,model,opts);
J = opDF3d(v,Q,model,opts);

x = randn(size(v)); x([1 end],:,:)= 0; x(:,[1 end],:) = 0; x(:,:,[1 end]) = 0; x = vec(x);
y = randn(size(D)) + 1i*randn(size(D)); y = vec(y);

s = J*x; s = real(y'*s);
t = J'*y; t = x'*t;
t = conj(t);
disp('Jacobian');
disp('<H*x,y>');
disp(num2str(s,'%3.15e'));
disp('<x,H''*y>');
disp(num2str(conj(t),'%3.15e'));
disp('Relative difference');
rel_diff = abs(s-t)/abs(s);
disp(rel_diff);


x = randn(size(v)); x([1 end],:,:)= 0; x(:,[1 end],:) = 0; x(:,:,[1 end]) = 0; x = vec(x);
y = randn(size(v)); y([1 end],:,:)= 0; y(:,[1 end],:) = 0; y(:,:,[1 end]) = 0; y = vec(y);

A = opH3d(v,Q,D,model,opts);
s = A*x; s = y'*s;
t = A'*y; t = x'*t;
disp('Hessian');
disp('<H*x,y>');
disp(num2str(s,'%3.15e'));
disp('<x,H''*y>');
disp(num2str(conj(t),'%3.15e'));
disp('Relative difference');
rel_diff = abs(s-t)/max(abs(s),abs(t));
disp(rel_diff);

