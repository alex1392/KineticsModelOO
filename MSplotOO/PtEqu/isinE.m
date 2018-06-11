function chk = isinE(E,obj)
if isVec(obj)
  for i = 1:size(obj,1)
    pt = obj(i,:);
    chk(i) = eqi(dot(E(1:3),pt),E(4));
  end
else
  for i = 1:numel(obj)
    if isa(obj,'SegL')
      sL = obj(i);
      chk(i) = all(isinE(E,sL.p));
    elseif isa(obj,'Line')
      L = obj(i);
      chk(i) = isperpend(L.vec, E(1:3)) && isinE(E,L.p);
    elseif isa(obj,'Face')
      F = obj(i);
      chk(i) = all(isinE(E,F.p));
    end
  end
end