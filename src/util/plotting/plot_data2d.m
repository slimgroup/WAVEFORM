function plot_data2d(D,model,opts)
% Plot velocity model
% vel_plot(v,model,opts)
%
% opts.
%    cmap      - colormap (default: 'default')
%    cbar      - show colorbar (default: false)
%    title_str - title string (default: [])
%    cax       - reference color axis (default: [])
%    slice_dir - for 3D models, direction to slice, one of 'x','y','z'
%    slice_depth - for 3D models, depth of slice
%
    
    if isfield(opts,'cmap'), cmap = opts.cmap; else cmap = seismic_colormap(200); end
    if isfield(opts,'cbar'), cbar = opts.cbar; else cbar = false; end
    if isfield(opts,'title_str'),title_str = opts.title_str; else title_str = ''; end
    if isfield(opts,'cax'), cax = opts.cax; else cax = []; end
    if isfield(opts,'new_fig'),new_fig = opts.new_fig; else new_fig = false; end
    if isfield(opts,'func'),func = opts.func; else func = @(x) x; end
    npts = 5;
    xlbl = 'source [km]'; ylbl = 'receiver [km]'; 
    xtick = floor(linspace(1,length(model.xsrc),npts)); 
    xticklabel = cellfun(@num2str,num2cell(model.xsrc(xtick)/1e3),'UniformOutput',false); 
    if max(xtick) > size(D,2),xtick(end) = size(D,2); end
    ytick = floor(linspace(1,length(model.xrec),npts)); 
    yticklabel = cellfun(@num2str,num2cell(model.xrec(ytick)/1e3),'UniformOutput',false);
    if max(ytick) > size(D,1),ytick(end) = size(D,1); end
    Ds = func(D);
    if new_fig, figure; end
    imagesc(Ds);
    colormap(cmap);
    title(title_str);
    if ~isempty(cax)
        caxis(cax);
    end
    xlabel(xlbl);    
    ylabel(ylbl);
    ax = gca;
    ax.XTick = xtick;
    ax.XTickLabel = xticklabel;
    ax.YTick = ytick;
    ax.YTickLabel = yticklabel;
    ax.FontSize = 14;    
    set(gcf,'Units','normalized');
    %set(gca,'LooseInset',get(gca,'TightInset'))
    %set(gcf,'OuterPosition',[0.1 0.1 0.5 0.5]);
    set(gcf,'PaperPositionMode','auto');
    if cbar, colorbar; end
end