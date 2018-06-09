function chk = isparallel(vA,vB,tol)
if nargin < 3
    tol = Const.FerroConst.Tol;
end
vA = vA/norm(vA);%unit vector
vB = vB/norm(vB);
chk = ~anyi(cross(vA,vB),tol);
end