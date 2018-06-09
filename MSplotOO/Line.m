classdef Line 
    
    properties (Hidden, Constant)
        Tol = Const.FerroConst.Tol;
    end
    
    properties %main
        p
        vec
    end
    
    properties (Hidden)
        long = 1;
    end
    
    methods
        %% constructor
        function L = Line(p, vec, varargin)
            if ~nargin
                return
            end
            L.p = p;
            L.vec = vec;
            if nargin > 2
                L = varassign(L, varargin);
            end
        end
        %% operational fcns
        function [rad,ang] = LvL(LA,LB)
          [rad,ang] = VvV(LA.vec,LB.vec);
        end
        
        function d = L_p(L,ps)
          for i = 1:size(ps,1)
            p = p(i,:);
            vA = p - L.p;
            vB = L.vec;
            d(i) = abs(cross(vA,vB))/norm(vB);
          end
        end
        
        function p = LxL(LAs,LBs)
            p = []; %only consider Point output
            nA = numel(LAs);
            nB = numel(LBs);
            for i = 1:nA
                for j = 1:nB
                    LA = LAs(i);
                    LB = LBs(j);
                    if isparallel(LA,LB)
                        p = []; 
                    else
                        A = [LA.vec; -LB.vec];
                        C = LB.p(1,:)-LA.p(1,:);
                        ts = C/A;
                        pA = LA.Lp(ts(1));
                        pB = LB.Lp(ts(2));
                        if eqi(pA,pB)
                            p = [p; pA];
                        end
                    end
                end
            end
            p = uniquei(p);
        end
        
        function p = LxE(L,E)
            p = []; %only consider Point output
            nL = numel(L);
            nE = size(E,1);
            for i = 1:nL
                for j = 1:nE
                    if isperpend(L(i).vec,E(j,1:3)) && isinE(E(j,:),L(i).p) %L in Equ
                        p = [p; L(i).p];
                    elseif ~isperpend(L(i).vec,E(j,1:3))
                        t = (E(j,4) - dot(L(i).p,E(j,1:3))) / dot(L(i).vec, E(j,1:3));
                        pt = L(i).Lp(t);
                        if L(i).isin(pt)
                            p = [p; pt];
                        end
                    end
                end
            end
            p = uniquei(p);
        end
        
        function h = draw(L,varargin)
            h = [];
            for i = 1:numel(L)
                u = L(i).vec / norm(L(i).vec);
                pA = L(i).p - 0.5*u*L(i).long;
                pB = L(i).p + 0.5*u*L(i).long;
                p = [pA;pB];
                h = [h, varassign(line(p(:,1),p(:,2),p(:,3)),varargin)];
            end
            axis(h(1).Parent,'equal');
        end
        
        function p = Lp(L,t)
            p = L.p(1,:)+L.vec*t;
        end
        %% is* methods
        function c = isin(L,p)
            n = size(p,1);
            c = false(1,n);
            for i = 1:n
              if eqi(p(i,:),L.p)
                c(i) = true;
              else
                v = p(i,:) - L.p;
                c(i) = isparallel(v,L.vec);
              end
            end
        end
        
        function c = isparallel(LA,LB)
            c = isparallel(LA.vec, LB.vec);
        end
        %% set, get methods
        function L = set.p(L,p)
            assert( isVec(p) || isempty(p), 'Invalid input "p".' ); 
            L.p = p;
        end
        
        function L = set.vec(L,vec)
            assert( isVec(vec) || isempty(vec), 'Invalid input "vec".' );
            L.vec = vec;
        end
    end
    
    methods (Static)
        function L = example
            L = Line([0 0 0],[1 1 1]);
        end
    end
end