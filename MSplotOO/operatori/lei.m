function chk = lei(A, B, Tol)
if nargin < 3
    Tol = Const.Tol;
end
chk = ~gti(A, B, Tol);
end