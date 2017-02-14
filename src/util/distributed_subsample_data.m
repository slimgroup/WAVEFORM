function y = distributed_subsample_data(x,iS,iF)
    nsrc = size(x,2); nfreq = size(x,3);
    [IS,IF] = ndgrid(iS,iF);
    I = sub2ind([nsrc,nfreq],vec(IS),vec(IF));
    y = x(:,I);
    spmd
        y = redistribute(y,codistributor('1d',2));
    end
    y = vec(y);
end