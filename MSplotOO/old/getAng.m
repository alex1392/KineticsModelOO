function ang = getAng(p, Equ)
% make vectors
n = size(p,1);
vec = [];
for i = 2:n
    vec = [vec; (p(i,:)-p(1,:))/ norm(p(i,:)-p(1,:))]; %unit vector
end
% judge orientation
ang = -Inf(n,1);
ang(2) = 0;
for i = 3:n
    c = cross(vec(1,:), vec(i-1,:));
    d = dot(vec(1,:), vec(i-1,:));
    ang(i) = acosd(d);
    if dot(c,Equ(1:3)) < 0 %right hand rule
        ang(i) = -ang(i);
    end
end
