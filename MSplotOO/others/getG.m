function G = getG(VarVol,EF,stress)
EP = getEProp(VarVol);
epsL = EP.s*stress + EP.d'*EF;
epsT = epsL + EP.epsR';
DT = EP.e*epsL + EP.k*EF + EP.PR';
U = 0.5*( epsL'*EP.c*epsL + EF'*EP.k*EF );
Omega =  - ( EF'*DT + stress'*epsT );
G = U + Omega;
end