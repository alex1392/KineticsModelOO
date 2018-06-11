function x = genvar(r)
  x = [];
  for a=1:6
    for b=1:6
      if r == 1
        x(end+1,:)=[a b];
      elseif r == 2
        for c=1:6
          for d=1:6
            x(end+1,:)=[a b c d];
          end
        end
      end
    end
  end
  x(1,:) = [];
  for i = length(x):-1:1
    a = x(i,:);
    if length(unique(a)) == 1
      x(i,:) = [];
    else
      for j = 1:i-1
        b = x(j,:);
        if r == 1
          b1 = [b(2) b(1)];
          b2 = [];
          b3 = [];
        elseif r == 2
          b1 = [b(2) b(1) b(4) b(3)];
          b2 = [b(3) b(4) b(1) b(2)];
          b3 = [b(4) b(3) b(2) b(1)];
        end
        if length(intersect(a,b)) == length(a) &&...
            isequal(a,b1) || isequal(a,b2) || isequal(a,b3)
          x(i,:) = [];
          break
        end
      end
    end
  end
end