function ary = insert(ary,idx,data)
  if isrow(ary)
    ary = [ary(1:idx-1), data, ary(idx:end)];
  else
    ary = [ary(1:idx-1,:); data; ary(idx:end,:)];
  end
  