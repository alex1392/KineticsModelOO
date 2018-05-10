function chk = gei(A, B, Tol)
if nargin < 3
    Tol = Const.Tol;
end
chk = ~ltii(A, B, Tol);
end