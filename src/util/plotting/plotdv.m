function func = plotdv(z,x,mref)
    
nz = length(z); maxz = max(z); nx = length(x); maxx = max(x);
sz = [nz,nx];

tov = @(x) x;
vec_to_str_cell = @(x) cellfun( @(s) num2str(s), mat2cell(vec(x)',1,ones(1,numel(x))),'UniformOutput',false);
plotopts = struct;
plotopts.xlabel = 'x [km] '; plotopts.ylabel = 'z [km] ';
plotopts.cmap = 'gray';

np = 5;

plotopts.ytick =  round(linspace(1,nz-1,np))+1;
zgrid = linspace(0,maxz,length(plotopts.ytick));
plotopts.yticklabel = vec_to_str_cell(round(zgrid/1000,1));

plotopts.xtick =  round(linspace(1,nx-1,np))+1;
xgrid = linspace(0,maxx,length(plotopts.xtick));
plotopts.xticklabel = vec_to_str_cell(round(xgrid/1000,1));
plotopts.axes_fontsize = 16;

func = @(x) multi_imagesc(plotopts,reshape(x,sz));


func(mref); ax = caxis;
plotopts.caxis = ax;
func = @(x) multi_imagesc(plotopts,reshape(tov(x),sz));