function out = is_var_local( Q )
%IS_VAR_LOCAL Returns true if variable is not distributed (or PCT is not
%installed)
    if license('test','parallel_toolbox')
        out = ~isdistributed(Q) && ~iscodistributed(Q);
    else
        out = true;
    end

end

