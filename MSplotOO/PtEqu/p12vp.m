function [rad,ang] = p12vp(p,E)
n = size(p,1);
rad = -Inf(n,1);
rad(2) = 0;
for i = 3:n
  vA = p(2,:)-p(1,:);
  vB = p(i,:)-p(1,:);
  rad(i) = VvV(vA,vB);
  if dot(cross(vA,vB),E(1:3)) < 0
    rad(i) = -rad(i);
  end
end
ang = rad*180/pi;