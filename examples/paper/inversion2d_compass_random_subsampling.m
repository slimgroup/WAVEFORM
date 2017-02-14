% This script performs stochastic FWI with box constraints, reproducing 
% the figures in Section 5.5.1 of 'A Unified 2D/3D Large Scale Software Environment for Nonlinear Inverse Problems' C. Da Silva and F. Herrmann (2017)

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

% Number of outer iterations (i.e., number of random redraws)
nouter_iter = 5;

% Number of inner iterations (i.e., number of effective passes through the data for each random draw of sources)
ninner_iter = 2;

% Size of frequency batch overlap
overlap_size = 2;

% Percentage of full sources to use initially
initial_src_size = 0.25;

% Percentage of full sources to use for an increment
increment_src_size = 0.2;


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
params.batch_mode = true;
opts = struct; 
opts.maxIter = nouter_iter;
opts.innerit = ninner_iter;

size_freq_batch = parpool_size(); overlap = overlap_size;
params.pdefunopts.src_est_mode = PDEopts.SRC_EST_NONE;
freq_partition = partition(nfreq,size_freq_batch,overlap);
part = splice(1:nsrc,parpool_size());
opts.bmax = vec(cellfun(@(x)length(x),part))';

opts.b0 = round(initial_src_size*opts.bmax);
opts.binc = round(increment_src_size*opts.bmax);
mlo = min(vec(m)); mhi = max(vec(m));

%% Vanilla FWI
model_err = [];
opts.true_sol = m;
tic;
mest = m0;
for j=1:size(freq_partition,1)
    fbatch = freq_partition(j,:);
    % Select only sources at this frequency batch
    srcfreqmask = false(nsrc,nfreq);
    srcfreqmask(:,fbatch) = true;
    params.srcfreqmask = srcfreqmask;
    obj = misfit_setup(mest,Q,Dobs,model,params);
    [mest,err] = minfunc_frugal(obj,mest,mlo,mhi,opts);
    model_err = [model_err err];
end
T = toc; disp(['Time elapsed ' num2str(T) 's']);
model_err = model_err(find(model_err));

%% Plot results
plot_v(m); 
plot_v(m0);
plot_v(mest);
figure; plot(model_err/norm(m),'LineWidth',3);
ylabel('Relative model error');
xlabel('Number of subproblems solved');
set(gca,'fontsize',16);
