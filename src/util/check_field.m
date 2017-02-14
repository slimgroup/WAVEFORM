% Checks if a specific string is a field of a given structure; if it is, returns
% its value, otherwise returns the default value.
%
% Use:
%   out = check_field(struct,str_field,d_val)
%
% Input:
%  struct    - structure or object to be checked
%  str_field - name of the field we are looking for; this must be a string.
%  d_val     - default value (in case the field is not found).
%
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
%-----------------------------------------------------------------------------
function out = check_field(s,str_field,d_val)

if ~isempty(s)
    if isa(s,'struct')
        if isfield(s,str_field)
            out = s.(str_field);
            return;
        end
    elseif isprop(s,str_field)
        out = s.(str_field); return;
    end
    
end
out = d_val;

end
