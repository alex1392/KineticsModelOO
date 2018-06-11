function C = unioni(A, B, dim, Tol)
    if nargin < 3
        dim = 1;
    end
    if nargin < 4
        Tol = Const.Tol;
    end
    C = uniquei([A; B], dim, Tol);
end