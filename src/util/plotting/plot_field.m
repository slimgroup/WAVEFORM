function func = plot_field(z,x,field,cax,zdir,xdir)
    
    if exist('cax','var')==0,cax = []; end 
    if exist('zdir','var')==0,zdir = 'z'; end
    if exist('xdir','var')==0,xdir = 'x'; end
    nz = length(z); maxz = max(z); nx = length(x); maxx = max(x);
    sz = [nz,nx];
    
    vec_to_str_cell = @(x) cellfun( @(s) num2str(s), mat2cell(vec(x)',1,ones(1,numel(x))),'UniformOutput',false);
    plotopts = struct;
    plotopts.axes_fontsize = 16;
    plotopts.xlabel = [xdir ' [m] ']; plotopts.ylabel = [zdir ' [m] '];

    np = 5;
    plotopts.ytick =  round(linspace(1,nz-1,np))+1;
    zgrid = linspace(0,maxz,length(plotopts.ytick));
    plotopts.yticklabel = vec_to_str_cell(floor(zgrid/100)*100);
    
    plotopts.xtick =  round(linspace(1,nx-1,np))+1;
    xgrid = linspace(0,maxx,length(plotopts.xtick));
    plotopts.xticklabel = vec_to_str_cell(floor(xgrid/100)*100);
    func = @(x) multi_imagesc(plotopts,reshape(x,sz));
    
    func(field); 
    if isempty(cax)
        cax = caxis;
    else
        caxis(cax);
    end
    plotopts.caxis = cax;
    func = @(x) multi_imagesc(plotopts,reshape(x,sz));