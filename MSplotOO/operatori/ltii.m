function chk = ltii(A, B, Tol)
if nargin < 3
    Tol = Const.Tol;
end
chk = all(A - B < -Tol);
end