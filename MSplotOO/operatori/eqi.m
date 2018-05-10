function chk = eqi(A, B, Tol)
if nargin < 3
    Tol = Const.Tol;
end
chk = all(abs(A-B) < Tol);
end