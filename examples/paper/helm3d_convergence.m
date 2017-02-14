%% Model setup 
% Velocity
v0 = 2000;

% Size of domain
L = [5000 5000 5000];

results_dir = '~/scratch/pdesoftware/helm_convergence/';

%% Set up source/receiver grid
% options for solving the helmholtz equations
lsopts = LinSolveOpts();
lsopts.tol = 1e-6;
lsopts.maxit = 10000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

model = struct;
model.o = [0 0 0];
model.xsrc = round(L(1)/2);
model.ysrc = round(L(2)/2);
model.zsrc = round(L(3)/2);

% source wavelet time shift
model.t0   = 0; 
model.unit = 'm/s';

% ricker wavelet peak frequency
model.f0 = 10;

%% 
num_wavelengths = [5 10 25 40 50]; nlambda = length(num_wavelengths);
npts_per_wavelength = [6 8 10]; 
np = length(npts_per_wavelength);

time = zeros(nlambda,np);
niter = zeros(nlambda,np);
npts_1d = zeros(nlambda,np);
opts.k = [3 5 3 5];
for i=1:nlambda
    for j=1:np
        freq = num_wavelengths(i)/max(L)*v0;
        d = v0/(npts_per_wavelength(j)*freq);
                        
        n = floor(L/d)+1;
        opts.dt = d*ones(1,3);
        opts.pml_max = inf;
        opts.pml = floor(v0/(freq*d));
        opts.solve_opts = lsopts;
        opts.mat_free = true;
        opts.n_threads = 10;
        opts.disp_output = true;
        opts.scheme = PDEopts.HELM3D_OPERTO27;
        v = v0*ones(n);
        model.d = d; model.n = n;
        
        [H,comp_grid] = discrete_helmholtz(v,model,freq,opts);        
        npts_1d(i,j) = comp_grid.nt(1);
        [xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
        sinc_window = 4;
        Ps_x = opSincInterp(model.xsrc,xt,sinc_window);
        Ps_y = opSincInterp(model.ysrc,yt,sinc_window);
        Ps_z = opSincInterp(model.zsrc,zt,sinc_window);
        Ps = opKron(Ps_z,Ps_y,Ps_x);
        
        q = Ps*( 1 );
        lsopts1 = copy(lsopts);
        if ~isprop(lsopts1,'output_freq')
            lsopts1.addprop('output_freq'); 
        end    
        lsopts1.output_freq = 1;
        lsopts1.precond = H.solve_opts.precond;
        
        tic,[~,res] = linearsolve(H,q,zeros(size(q)),1,lsopts1); 
        T = toc;
        disp(['nlam ' num2str(num_wavelengths(i)) ' ppw ' num2str(npts_per_wavelength(j))  ' niter ' num2str(length(res)+1) ' time : ' num2str(T) 's']);
        time(i,j) = T;
        niter(i,j) = length(res)+1;
        
    end
end
