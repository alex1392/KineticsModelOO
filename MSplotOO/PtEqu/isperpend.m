function chk = isperpend(vA,vB,tol)
if nargin < 3
    tol = Const.Tol;
end
vA = vA/norm(vA);%unit vector
vB = vB/norm(vB);
chk = eqi(dot(vA,vB),0,tol);
end