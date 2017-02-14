function func = plotv(z,x,mref,unit)
    
nz = length(z); maxz = max(z); nx = length(x); maxx = max(x);
sz = [nz,nx];
if exist('unit','var')==0, unit = 's2/km2'; end
if strcmp(unit,'s2/km2')
    tov = @(x) (1e3*(x).^(-1/2));
else
    tov = @(x) x;
end

vec_to_str_cell = @(x) cellfun( @(s) num2str(s), mat2cell(vec(x)',1,ones(1,numel(x))),'UniformOutput',false);
plotopts = struct;
plotopts.xlabel = 'x [km] '; plotopts.ylabel = 'z [km] ';

np = 5;
plotopts.ytick =  round(linspace(1,nz-1,np))+1;
zgrid = linspace(0,maxz,length(plotopts.ytick));
plotopts.yticklabel = vec_to_str_cell(round(zgrid/1000,1));

plotopts.xtick =  round(linspace(1,nx-1,np))+1;
xgrid = linspace(0,maxx,length(plotopts.xtick));
plotopts.xticklabel = vec_to_str_cell(round(xgrid/1000,1));
plotopts.axes_fontsize = 16;
func = @(x) multi_imagesc(plotopts,reshape(tov(abs(x)),sz));

func(mref); ax = caxis;
plotopts.caxis = ax;

func = @(x) multi_imagesc(plotopts,reshape(tov(abs(x)),sz));