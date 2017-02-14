% This script performs straightforward, frequency continuation FWI, reproducing 
% the figures in Section 5.2 of 'A Unified 2D/3D Large Scale Software Environment for Nonlinear Inverse Problems' C. Da Silva and F. Herrmann (2017)

% Minimum frequency [Hz]
min_freq = 3;

% Maximum frequency [Hz]
max_freq = 18;

% Number of frequencies
nfreq = max_freq-min_freq+1;

% Number of sources
nsrcs = 100;

% Number of parallel workers
nworkers = 4;

% Number of smoothing of the true model
nsmooth = 100;

% Number of optimization iterations
niter = 20;

% Size of frequency batch overlap
overlap_size = 2;

%% Load model and generate data
params = struct; 
params.nsrcs = nsrcs; 
params.minFreq = min_freq; 
params.nFreq = nfreq;
params.maxFreq = max_freq;

model = loadModel('compass_small',params);
m = vec(model.m);
nsrc = length(model.xsrc); nrec = length(model.xrec); nfreq = length(model.freq); 
nz = model.n(1); nx = model.n(2);
Q = speye(nsrc);

if parpool_size()==0, pool = parpool(nworkers); end
Dobs = F(m,Q,model);
plot_v = plotv(model.z,model.x,m);

%% FWI setup
m0 = opKron(opSmooth(nx,nsmooth),opSmooth(nz,nsmooth))*m;
params = struct;
params = default_fwi_params2d(params);
params.wri = false;
opts = struct; 
opts.maxIter = niter;
size_freq_batch = parpool_size(); overlap = overlap_size;
params.pdefunopts.src_est_mode = PDEopts.SRC_EST_NONE;
freq_partition = partition(nfreq,size_freq_batch,overlap);
mlo = min(vec(m)); mhi = max(vec(m));

%% Vanilla FWI
tic,
mest = m0;
for j=1:size(freq_partition,1)
    fbatch = freq_partition(j,:);
    % Select only sources at this frequency batch
    srcfreqmask = false(nsrc,nfreq);
    srcfreqmask(:,fbatch) = true;
    params.srcfreqmask = srcfreqmask;
    obj = misfit_setup(mest,Q,Dobs,model,params);
    mest = minConf_TMP(obj,mest,mlo,mhi,opts);    
end
mest_vanilla = mest;
T = toc; disp(T);

%% Display results
plot_v(m); title('True model');
plot_v(m0); title('Initial model');
plot_v(mest); title('Estimated model');

