function L = ExE(EquA, EquB, tol)
if nargin < 3
    tol = Const.FerroConst.Tol;
end
A = [EquA(1:3); EquB(1:3)];
if rank(A, tol) < 2
    L = [];
else
    B = [ EquA(4); EquB(4) ];
    L = Line( (A\B)', cross( A(1,:), A(2,:) ) );
end