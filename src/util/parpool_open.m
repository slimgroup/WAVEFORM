function [ pool ] = parpool_open( nworkers )
%PARPOOL_OPEN Opens default parallel pool of given size
    pool = parpool(nworkers)

end

