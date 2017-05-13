function ps = parpool_size()
% parpool_size - return number of workers in parallel pool
% parpool_size returns 0 id parallel pool is not opened
    if license('test','parallel_toolbox')
        pool = gcp('nocreate');
        if isempty(pool)
            ps = 0;
        else
            ps = pool.NumWorkers;
        end
    else
        ps = 0;
    end
end
