function chk = anyi(obj, Tol)
if nargin < 2
    Tol = Const.Tol;
end
chk = any(abs(obj)>Tol);