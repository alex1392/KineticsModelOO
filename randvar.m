function [vars,dofs] = randvar(n)
vars = {};
for i = 1:n
  first = randi(6);
  var = first;
  switch randi(8)
    case 1 %1112
      var(2:3) = first;
      var(4) = get180(first);
    case 2 %1213
      var(2) = get180(first);
      var(3) = first;
      tmp = setdiff(1:6,var(1:2));
      odr = randperm(length(tmp));
      var(4) = tmp(odr(1));
    case 3 %1221
      var(4) = first;
      var(2:3) = get180(first);
    case 4 %1325
      var(3) = get180(first);
      tmp1 = setdiff(1:2:5,var([1 3]));
      odr = randperm(length(tmp1));
      tmp1 = tmp1(odr);
      tmp2 = randi(1,1,2);
      tmp = tmp1+tmp2;
      var(2) = tmp(1);
      var(4) = tmp(2);
    case 5 %1234
      var(2) = get180(first);
      tmp = setdiff(1:2:5,var(1:2));
      odr = randperm(length(tmp));
      var(3) = tmp(odr(1));
      var(4) = get180(var(3));
    case 6 %1314
      var(3) = first;
      tmp = setdiff(1:2:5,[first get180(first)]);
      odr = randperm(length(tmp));
      var(2) = tmp(odr(1));
      var(4) = get180(var(2));
    case 7 %1324
      var(3) = get180(first);
      tmp = setdiff(1:2:5,var([1 3]));
      odr = randperm(length(tmp));
      var(2) = tmp(odr(1));
      var(4) = get180(var(2));
    case 8 %1342
      var(4) = get180(first);
      tmp = setdiff(1:2:5,var([1 4]));
      odr = randperm(length(tmp));
      var(2) = tmp(odr(1));
      var(3) = get180(var(2));
  end
  vars{end+1} = var;
end
dofs = repmat({[0.5 0.5]},[1 n]);


function m = get180(n)
switch n
  case 1
    m = 2;
  case 2
    m = 1;
  case 3
    m = 4;
  case 4
    m = 3;
  case 5
    m = 6;
  case 6
    m = 5;
end
    
    
    