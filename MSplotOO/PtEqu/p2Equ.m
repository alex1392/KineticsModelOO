function Equ = p2Equ(p)
vec = [p(1,:) - p(2,:); p(1,:) - p(3,:)];
eqN = cross(vec(1,:), vec(2,:));
d = dot(eqN, p(1,:));
Equ = [eqN,d]/norm(eqN);
end