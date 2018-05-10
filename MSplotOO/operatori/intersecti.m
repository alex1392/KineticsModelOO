function [C, ia, ib] = intersecti(A, B, dim, Tol)
    if nargin < 3
        dim = 1;
    end
    if nargin < 4
        Tol = Const.Tol;
    end
    if dim == 2
        A = A';
        B = B';
    end
    if isa(A,'Point') || isa(A,'Plane')
        nA = A.n;
    else
        nA = size(A,1);
    end
    if isa(B,'Point') || isa(B,'Plane')
        nB = B.n;
    else
        nB = size(B,1);
    end
    C = [];  ia = [];  ib = [];
    for ra = 1:nA
        for rb = 1:nB
            if eqi(A(ra,:), B(rb,:), Tol)
                C = [C; A(ra,:)];
                ia = [ia; ra];
                ib = [ib; rb];
            end
        end
    end
    if dim == 2
        C = C';  ia = ia';  ib = ib';
    end
end