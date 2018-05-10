function c = isEqu( Equ, n )
    %%
    if nargin < 2
        n = size(Equ,1);
    end
    if isnumeric(Equ) && isequal( size(Equ), [n,4] )
        c = 1;
    else
        c = 0;
    end
end