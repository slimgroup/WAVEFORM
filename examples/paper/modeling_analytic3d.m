%% Model geometry
n1d = 201;
nrec = 101;
model.o = [0 0 0];
model.n = [n1d n1d n1d];
model.d = 2000./ (model.n-1);
[x,y,z] = odn2grid(model.o,model.d,model.n);

freq = 8;
model.freq = freq;
model.f0 = 10;
model.t0 = 0;
model.unit = 'm/s';
model.zsrc = 20;
model.xsrc = 0;
model.ysrc = 100;
model.zrec = 40;
model.xrec = linspace(0,max(x),nrec);
model.yrec = 0;
Q = speye(length(model.xsrc)*length(model.ysrc));

%% Velocity + wavefields
v0 = 2000;
vc = v0*ones(model.n);
lsopts = LinSolveOpts();
lsopts.tol = 1e-10;
lsopts.maxit = 2000;
lsopts.maxinnerit = 5;
lsopts.solver = LinSolveOpts.SOLVE_FGMRES;
lsopts.precond = LinSolveOpts.PREC_MLGMRES;

opts = struct;
opts.free_surface = false;
opts.disp_output = true;
opts.scheme = PDEopts.HELM3D_OPERTO27;
opts.solve_opts = lsopts;
opts.pml = 30;

[Hk,comp_grid] = discrete_helmholtz(vc,model,freq,opts);

[xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
[x,y,z] = odn2grid(model.o,model.d,model.n);

Ptr = opKron(opExtension(model.n(3),[opts.pml opts.pml],0),opExtension(model.n(2),[opts.pml opts.pml],0),opExtension(model.n(1),[opts.pml opts.pml],0));

Ps = opInterp('sinc',model.xsrc,xt,model.ysrc,yt,model.zsrc,zt);

q = Ps*Q(:,1);

tic,uc = Hk\q;T3d = toc;disp(T3d);
%close all;

ixsrc = model.xsrc/model.d(1) + 1;
iysrc = model.ysrc/model.d(2) + 1;
izsrc = model.zsrc/model.d(3) + 1;

u = reshape(Ptr'*uc,model.n);

ax = [-1 1]*0.25;

plotRx = plot_field(y,z,real(u(ixsrc,:,:)),ax,'y','z'); title('Computed solution : real');
plotIx = plot_field(y,z,imag(u(ixsrc,:,:)),ax,'y','z'); title('Computed solution : imag');

    
%% Green's function
[~,I] = max(abs(vec(q))); 
[ix,iy,iz] = ind2sub(comp_grid.nt,I);
[XT,YT,ZT] = ndgrid(xt,yt,zt);
R = ((XT-xt(ix)).^2 +(YT-yt(iy)).^2 + (ZT-zt(iz)).^2 ).^(1/2);
w = (2*pi*freq/v0);
G = prod(model.d)*exp(1i*w*R)./(4*pi*R);

G = reshape(Ptr'*vec(G),model.n);
plotRx(real(G(ixsrc,:,:))); title('Greens function : real');

plotRx(real(G(ixsrc,:,:)-u(ixsrc,:,:))); title('Analytic - computed : real, x10');
caxis(ax*1e-1);

plotIx(imag(G(ixsrc,:,:))); title('Greens function : imag');
  
plotIx(imag(G(ixsrc,:,:)-u(ixsrc,:,:))); title('Analytic - computed : imag, x10');
caxis(ax*1e-1);