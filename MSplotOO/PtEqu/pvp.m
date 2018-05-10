function [rad,ang] = pvp(p,pA,pB)
vA = (pA - p)/norm(pA - p);
vB = (pB - p)/norm(pB - p);
[rad,ang] = VvV(vA,vB);