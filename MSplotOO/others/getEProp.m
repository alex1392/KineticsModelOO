function EffProp = getEProp(VarVol)
EffProp = EProp;
Vol = sum(VarVol);
for V = 1:Const.vn
  f = VarVol(V) / Vol;
  EffProp = EffProp + Const.MProp(V)*f;
end