% This script performs a simple conductivity inversion for a square ananomaly on a constant background, reproducing the figures in Section 5.4 of 'A Unified 2D/3D Large Scale Software Environment for Nonlinear Inverse Problems' C. Da Silva and F. Herrmann (2017)

% Number of points in each direction
n = [201,201];

% Grid spacing
d = 25;

% Number of points to restrict to a constant value at the top + bottom of the model
np = 40;

% Permittivity of free space
eps0 = 8.854*1e-12;

% Size of square anomaly
s = 20;

epsilon = ones(n);
epsilon(round(n(1)/2)-s:round(n(1)/2)+s,round(n(2)/2)-s:round(n(2)/2)+s) = 1.2;

%% Construct model + generate data
model = struct;
model.n = n;
model.o = [0 0];
model.d = [d d];
model.zsrc = np*model.d(1);
model.xsrc = linspace(0,model.d(2)*model.n(2),100);
model.zrec = model.d(1)*(model.n(1)-np);
model.xrec = linspace(0,model.d(2)*model.n(2),model.n(2));
model.freq = 1;
model.f0 = 10;
model.t0 = 0;
nsrc = length(model.xsrc); nrec = length(model.xrec);
params = struct;
params.wri = false;
params = default_fwi_params2d(params);
params.pdefunopts.helm_scheme = PDEopts.POISS2D_FV;
params.pdefunopts.helm_pml = 0;
params.hessian = PDEopts.HESS_GN;

Q = -speye(nsrc)/eps0;

D = reshape(F(epsilon,Q,model,params),nsrc,nrec);

epsilon0 = ones(n);
D0 = reshape(F(epsilon0,Q,model,params),nsrc,nrec);

%% Gauss-Newton method
elo = ones(n);
ehi = 2*max(vec(epsilon))*ones(n);
ehi(1:np+5,:) = 1; ehi(end-(np+5):end,:) = 1;    

opts = struct;
opts.maxIter = 5;

obj = misfit_setup(epsilon0,Q,D,model,params);
e = vec(epsilon0);
for i=1:5
    [f,g,h] = obj(e);
    disp(f);
    q = @(x) quad_model(g,h,x,e);
    e = minConf_TMP(q,vec(e),vec(elo),vec(ehi),opts);
end

%% Plot results
[z,x] = odn2grid(model.o,model.d,model.n);
p = plotv(z,x,epsilon,'m/s'); title('True model');
p(epsilon0); title('Initial model');
p(e); title('Estimated model');
