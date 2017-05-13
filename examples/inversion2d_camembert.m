%% Camembert example 
% This script is based on the famous `Camembert' examples from 
% O. Gauthier, J. Virieux, and A. Tarantola. Two-dimensional nonlinear
% inversion of seismic waveforms: Numerical results. Geophysics 51, 1387-1403 (1986)
% A simple example illustrating the inversion framework on the famous
% Camembert model, with reflection and transmission source/receiver configurations

%% define velocity model

% grid
z = 0:10:1000;
x = 0:10:1000;
[o,d,n] = grid2odn(z,x);
[zz,xx] = ndgrid(z,x);

% background velocity 2500 [m/s]
v0 = 2500 + 0*xx;

% circular perturbation with radius 250 m and strenth 10\%
dv = 0*xx; dv((xx-500).^2 + (zz-500).^2 <= 250^2) = .1*2500;

% plot
figure;imagesc(x,z,v0 + dv);xlabel('x [m]');ylabel('z [m]');title('velocity [m/s]');colorbar;


%% setup modeling parameters for reflection experiment

% grid
model_r.o = o; 
model_r.d = d;
model_r.n = n;

% number of points for absorbing boundary
model_r.nb = [20 20];

% frequencies [5, 10, 15] Hz.
model_r.freq = [5 10 15]; nfreq = length(model_r.freq);

% Ricker wavelet with peak frequency f0 and phase shift t0
model_r.f0 = 10;
model_r.t0 = 0;

% source and receiver locations
model_r.zsrc = 10;
model_r.xsrc = 0:100:1000; nsrc = length(model_r.xsrc);
model_r.zrec = 10;
model_r.xrec = 0:10:1000;  nrec = length(model_r.xrec);

% source matrix, each column is a source function defined on the source
% grid [model.zsrc, model.xsrc].
Q = speye(nsrc);

% squared slowness [km^2/s^2]
m = 1e6./(v0(:) + dv(:)).^2;

%% create data
D = F(m,Q,model_r);

%% inversion
opts = struct;
opts.wri = false;

% initial model
m0 = 1e6./v0(:).^2;
obj = misfit_setup(m0,Q,D,model_r,opts);

optim_opts = struct;
optim_opts.maxIter = 20;

m_r = minConf_TMP(obj,m0,min(m),max(m),optim_opts);

% plot
figure;imagesc(x,z,reshape(1e3*(m_r.^(-1/2)),n));xlabel('x [m]');ylabel('z [m]');title('reconstruction from reflection data');

%% setup modeling parameters for transmission experiment

% grid
model_t.o = o; 
model_t.d = d;
model_t.n = n;

% number of points for absorbing boundary
model_t.nb = [20 20];

% frequencies [5, 10, 15] Hz.
model_t.freq = [5 10 15]; nfreq = length(model_t.freq);

% Ricker wavelet with peak frequency f0 and phase shift t0
model_t.f0 = 10;
model_t.t0 = 0;

% source and receiver locations
model_t.zsrc = 0:100:1000; nsrc = length(model_t.zsrc);
model_t.xsrc = 50;
model_t.zrec = 0:10:1000;  nrec = length(model_t.zrec);
model_t.xrec = 950;

% source matrix, each column is a source function defined on the source grid [model.zsrc, model.xsrc].
Q = speye(nsrc);

% squared slowness [km^2/s^2]
m = 1e6./(v0(:) + dv(:)).^2;

% create data
D_t = F(m,Q,model_t);

%% inversion

% initial model
m0 = 1e6./v0(:).^2;

obj = misfit_setup(m0,Q,D_t,model_t,opts);
m_t =  minConf_TMP(obj,m0,min(m),max(m),optim_opts);

% plot
figure;imagesc(x,z,reshape(1e3*(m_t.^(-1/2)),n));
xlabel('x [m]');ylabel('z [m]');title('reconstruction from transmission data');

