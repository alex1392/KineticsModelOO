function [A, ia] = setdiffi(A, B, dim, Tol)
    %% SETDIFF  Kill the same elements in objA compared to objB alone dim.
    if nargin < 3
        dim = 1;
    end
    if nargin < 4
        Tol = Const.FerroConst.Tol;
    end
    [~, ia, ~] = intersecti(A, B, dim, Tol);
    A(ia,:) = [];
end
    