function chk = alli(obj, Tol)
if nargin < 2
    Tol = Const.Tol;
end
chk = all(abs(obj)>Tol);