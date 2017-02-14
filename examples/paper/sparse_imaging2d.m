% This script performs sparse imaging with randomized subsampling using the 
% linearized Bregman algorithm, reproducing the figures in Section 5.3 of 'A Unified 2D/3D Large Scale Software Environment for Nonlinear Inverse Problems' C. Da Silva and F. Herrmann (2017) 
% 
% This example requires the installation of Curvelab (http://www.curvelet.org/) and LSMR (https://web.stanford.edu/group/SOL/software/lsmr/) 

% Random seed
seed = 128575234;

% Minimum frequency [Hz]
minFreq = 3; 

% Maximum frequency [Hz]
maxFreq = 12; 

% Number of frequencies
nFreq = 40;

% Number of sources
nsrc = 300;

% Amount of smoothing of the true model
nsmooth = 40;

% Number of passes through the total data
nDataPasses = 10;

% Fraction of sources to use for each random batch
psrc = 0.1;

% Number of parallel workers
nworkers = 10;

%% Load model and generate data
params = struct; 
params.nsrcs = nsrc; 
params.minFreq = minFreq; 
params.nFreq = nFreq;
params.maxFreq = maxFreq;
model = loadModel('compass_small',params);
m = vec(model.m);
nsrc = length(model.xsrc); nrec = length(model.xrec); nfreq = length(model.freq); 
nz = model.n(1); nx = model.n(2);
Q = speye(nsrc);
if parpool_size()==0, pool = parpool(nworkers); end
Dobs = F(m,Q,model);

m0 = opKron(opSmooth(nx,nsmooth),opSmooth(nz,nsmooth))*m;
params = struct;
params = default_fwi_params2d(params);
params.wri = false;

%% Construct operators 
dm = m-m0;
A_full = oppDF(m0,Q,model,params);
b_full = A_full*dm;
b_full = reshape(b_full,nrec,nsrc,nfreq);
NBscales = max(1,ceil(log2(min(model.n(1:2))) - 3));
NBangles = 16;
Finest = 1;
Ttype = 'ME';
IS_real = 1;
C = opCurvelet(model.n(1),model.n(2),NBscales,NBangles,Finest,Ttype,IS_real);

%% Linearized Bregman with randomized subsampling
softThreshold = @(x,lambda) max(abs(x)-lambda,0).*sign(x);
rng(seed);
x = 0*dm;
z = x;
t = psrc*parpool_size()/nfreq;
T = round(nDataPasses/t);
rand_subset = @(N,M) sort(randperm(N,M),'ascend');
res = zeros(T,1);
for i=1:T
    % Subsample A, b
    Is = rand_subset(nsrc,round(psrc*nsrc));
    If = rand_subset(nfreq,parpool_size());
    srcfreqmask = false(nsrc,nfreq);
    srcfreqmask(Is,If) = true;
    params.srcfreqmask = srcfreqmask;
    A = oppDF(m0,Q,model,params);
    b = distributed_subsample_data(b_full,Is,If);
    if i==1
        % Initial guess is zero
        r = -b;
    else
        r = A*x-b;
    end    
    ATr = A'*r;
    t = norm(r)^2/norm(ATr)^2;
    if i==1
        lambda = 0.1*max(abs(C*(t*ATr)));
    end
    z = z - t*ATr;
    x = C'*softThreshold(C*z,lambda);
    res(i) = norm(x-dm)/norm(dm);
    disp(['itr ' num2str(i) ' : ' num2str(res(i))]);
end

%% Least-squares inversion with the full data 
x_full = lsmr(A_full,vec(b_full),0,1e-6,[],[],nDataPasses,0,true);

%% Plot results
plot_v = plotdv(model.z,model.x,dm);
caxis([-1 1]*0.02); ax = caxis; title('True perturbation');

plot_v(x); caxis(ax); title('Linearized Bregman image');
plot_v(x_full); caxis(ax); title('Full data image');


