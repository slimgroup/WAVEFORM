function model = loadModel(modelName,params)
% loadModel - 
%
% Input parameters
%    .nsrcs 
%    .minFreq
%    .maxFreq
%    .nFreq
%    .z
%    .x 
%    .transmission

nsrcs = params.nsrcs; 
minFreq = params.minFreq; maxFreq = params.maxFreq; 
nFreq = params.nFreq;
if isfield(params,'transmission'),transmission = params.transmission; else transmission = false; end

switch modelName
  case {'slope','camembert','bump'}
    z = params.z; x = params.x; 
    nz = length(z); nx = length(x);
        
    % Model generation
    x2 = linspace(0,1,nx); z2 = linspace(0,1,nz);
    [Z,X] = ndgrid(z2,x2);
    m = 1e6/2500^2 * ones(nz,nx);
    if strcmp(modelName,'slope')
        m( 8*(Z-0.4) + X >= 0 ) = 1e6/2700^2;
        m( 8*(Z-0.8) + X >= 0 ) = 1e6/3000^2;
    elseif strcmp(modelName,'camembert')
       m( (Z- 0.5).^2 + (X-0.5).^2 <= 0.25^2 ) = 1e6/(1.1*2500)^2; 
    elseif strcmp(modelName,'bump')
        m = 2500 * ones(nz,nx);
        R = params.r;
        m = m + 1000*exp( -((Z-0.5).^2 + (X-0.5).^2)/R^2 );
        m( 8*(Z-0.85) + X >= 0 & 8*(Z-0.88) +X <= 0 ) = 3000;    
        m = 1e6 * (m.^(-2));
    end
    model.m = m;
    model.nz = nz; model.nx = nx;    
    
    % Recording paramters
    t    = 0:.004:2; 
    freq = 0:1/t(end):.5/(t(2)-t(1));
    [o,d,n] = grid2odn(z,x,0);
    model.o = o; model.d = d; model.n = n;
    model.nb = [20 20 0];
    [~,minIf] = min(abs(freq-minFreq)); [~,maxIf] = min(abs(freq-maxFreq));
    If = round(linspace(minIf,maxIf,nFreq));
    
    model.freq = freq(If);            
    if transmission
        offset = round(0.15*length(x));
        model.xsrc = x(offset);
        model.zsrc = linspace(0,max(z),nsrcs);
        model.xrec = x(end-offset);
        model.zrec = z;
    else
        model.zsrc = 10;
        model.xsrc = linspace(0,max(x),nsrcs);  
        model.zrec = 10;
        model.xrec = linspace(0,max(x),nx); 
    end
    model.f0   = 15;
    model.t0   = 0;
    
  case 'compass_small'
    modelFile = 'BG450_model.mat';
            
    load(modelFile);
    
    v = v*1e3;

    v = reshape(v,n);
    n = size(v);
    
    n = [n,1]; d = [d,0]; o = [o,0];
    [z,x] = odn2grid(o,d,n);
        
    
    model = struct;
    model.o = o;
    model.d = d;
    model.n = n;
    model.nb = [40 40 0];
    model.z = z; model.x = x;
    t    = 0:.004:2; nt = length(t);
    freq = 0:1/t(end):.5/(t(2)-t(1));
    [~,minIf] = min(abs(freq-minFreq)); [~,maxIf] = min(abs(freq-maxFreq));
    If = round(linspace(minIf,maxIf,nFreq));
    
    model.freq = freq(If);       
    minFreq = min(model.freq); maxFreq = max(model.freq);
    model.zsrc = 10;
    model.xsrc = linspace(0,max(x),nsrcs); 
    model.zrec = 10;
    model.xrec = linspace(0,max(x),n(2)); model.nrec  = length(model.xrec);
    model.f0   = 15;
    model.t0   = 0;
    model.m = 1e6./(v.^2);  
end
    
end