function chk = nei(A, B, Tol)
if nargin < 3
    Tol = Const.Tol;
end
chk = ~eqi(A, B, Tol);
end