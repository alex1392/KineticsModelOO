function [rad,ang] = VvV(vA,vB)
vA = vA/norm(vA);
vB = vB/norm(vB);
rad = acos(dot(vA,vB));
ang = rad*180/pi;