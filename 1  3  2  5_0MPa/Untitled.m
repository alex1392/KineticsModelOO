f = Face([1 0 0;0 1 0;0 0 1]);

on = MStruct.empty;
below = MStruct.empty;
for i = MST.findleaves
  if all(f.ison(MST.get(i).p))
    on(end+1) = MST.get(i);
  elseif ~any(f.ison(MST.get(i).p))
    below(end+1) = MST.get(i);
  else
    [Pon,Pbelow,~] = MST.get(i).Cut(f.Equ);
    on(end+1) = Pon;
    below(end+1) = Pbelow;
  end
end
on.draw;
figure
below.draw;