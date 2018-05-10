function [Lia, Locb] = ismemberi(A, B, dim, Tol)
if nargin < 3
    dim = 1;
end
if nargin < 4
    Tol = Const.Tol;
end
if dim == 2
    A = A';  B = B';
end
[~, ia, ib] = intersecti(A, B, dim, Tol);
Lia = false(size(A,1),1);
Lia(ia) = true;
Locb = int8(Lia);
Locb(ia) = ib;
if dim == 2
    Lia = Lia';  Locb = Locb';
end
end