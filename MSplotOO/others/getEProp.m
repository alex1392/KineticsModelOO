function EffProp = getEProp(VarVol)
EffProp = EProp;
Vol = sum(VarVol);
for V = 1:Const.FerroConst.vn
  f = VarVol(V) / Vol;
  EffProp = EffProp + Const.FerroConst.MProp(V)*f;
end