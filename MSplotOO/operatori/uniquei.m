function obj = uniquei(obj,dim,Tol)
if isa(obj,'Point')
  obj = obj.vec;
  ispt = 1;
else
  ispt = 0;
end
if nargin < 2
  dim = 1;
end
if nargin < 3
  Tol = Const.FerroConst.Tol;
end
if dim == 2
  obj = obj';
end
i = 1;
while i <= size(obj,1) %must check it's number every time
  for j = size(obj,1) : -1 : i+1
    if eqi(obj(i,:), obj(j,:), Tol)
      obj(j,:) = [];
    end
  end
  i = i + 1;
end
if dim == 2
  obj = obj';
end
if ispt
  obj = Point(obj);
end
end
