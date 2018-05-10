function chk = isonE(E,obj)
if isVec(obj)
    np = size(obj,1);
    chk = false(1,np);
    for i = 1:np
        p = obj(i,:);
        chk(i) = dot(p, E(1:3)) > E(4);
    end
elseif isa(obj, 'SegL') || isa(obj, 'Face')
    n = numel(obj);
    chk = false(1,n);
    for i = 1:n
        chk(i) = all(F.ison(obj(i).p));
    end
end
end