function obj = varassign(obj,var)
if nargin < 2 || isempty(var)
    return
end
nvarin = numel(var);
assert(~mod(nvarin,2),'Field and value input arguments must come in pairs.');
for i = 1:2:nvarin-1
    obj.(var{i}) = var{i+1};
end