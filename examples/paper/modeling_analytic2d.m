%% Demonstrates the accuracy of the 2D Helmholtz discretization for a constant velocity model
    
%% Model geometry
nz = 201; nx = 201;
nsrc = 50; nrec = 150;
model.o = [0 0 0];
model.n = [nz nx 1];
model.d = [10 10 1];
[z,x] = odn2grid(model.o,model.d,model.n);
dx = model.d(2);
model.freq = 4:8;
model.f0 = 10;
model.t0 = 0;
model.unit = 'm/s';
model.zsrc = 100;
model.xsrc = linspace(4*dx,max(x)-4*dx,nsrc);
model.zrec = 100;
model.xrec = linspace(4*dx,max(x)-4*dx,nrec);
model.nb = [20 20 0];

Q = speye(nsrc);
    
    
%% Velocity + wavefields
v0 = 2000;
vc = v0*ones(model.n);

Pext = opKron(opExtension(model.n(2),[model.nb(2) model.nb(2)]),opExtension(model.n(1),[model.nb(1) model.nb(1)]));
Ptr = opKron(opExtension(model.n(2),[model.nb(2) model.nb(2)],0),opExtension(model.n(1),[model.nb(1) model.nb(1)],0));

vx = Pext*vec(vc);
ot = model.o - model.nb.*model.d;
nt = model.n + 2*model.nb; 

k = 1;
Hk = Helm2D_opt(vx,model.d,nt,vec(model.nb(1:2))*ones(1,2),'m/s',model.freq(k),model.f0);
[z,x] = odn2grid(model.o,model.d,model.n);
[zt,xt] = odn2grid(ot,model.d,nt);
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc))';

q = Ps*Q(:,1);
uc = Hk\q;

uc = Ptr'*vec(uc);

cax = [-1 1]*10;

plotR = plot_field(z,x,real(uc),cax); title('Computed field : real');
plotI = plot_field(z,x,imag(uc),cax); title('Computed field : imag');

%% Green's function
[~,I] = max(abs(vec(q))); [iz,ix] = ind2sub(nt,I);
[ZT,XT] = ndgrid(zt,xt);
R = ((ZT-zt(iz)).^2 + ((XT-xt(ix)).^2)).^(1/2);
t = (2*pi*model.freq(k)/v0)*R;
G = -1i/4*besselh(0,t);
G(iz,ix) = -0.5-0.25*1i;
G = prod(model.d)* conj(G);

G = Ptr'*vec(G);
plotR(real(G)); title('Greens function : real');
plotR(real(G-uc)); title('Analytic - computed : real, x100');
caxis(cax*1e-2); 

plotI(imag(G)); title('Greens function : imag');
plotI(imag(G-uc)); title('Analytic - computed : imag, x100');
caxis(cax*1e-2);



