function chk = isVec(vec, n)
    if nargin < 2
        n = size(vec,1);
    end
    chk = isnumeric(vec) && isequal( size(vec), [n,3] );
end