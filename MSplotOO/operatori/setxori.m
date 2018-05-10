function [C, xa, xb] = setxori(A, B, dim, Tol)
    %% SETXOR  Find the elements in objA and objB that are not in their intersection along dim.
    if nargin < 3
        dim = 1;
    end
    if nargin < 4
        Tol = Const.Tol;
    end
    xa = 1:size(A,1);
    xb = 1:size(B,1);
    [~, ia, ib] = intersecti(A, B, dim, Tol);
    A(ia,:) = [];  xa(ia) = [];
    B(ib,:) = [];  xb(ib) = [];
    C = [A; B];
end