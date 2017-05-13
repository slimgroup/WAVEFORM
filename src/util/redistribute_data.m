function Dobs = redistribute_data( Dobs )
spmd,
    Dobs = redistribute(Dobs,codistributor1d(2,codistributor1d.unsetPartition,size(Dobs)));
end

end

