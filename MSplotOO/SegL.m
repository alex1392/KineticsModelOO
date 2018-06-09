classdef SegL
    properties (Hidden, Constant)
        Tol = Const.FerroConst.Tol;
    end
    
    properties
        %% main porperties
        p
        vec
        long
    end
    
    methods
        %% constructor
        function sL = SegL(p, varargin)
            if ~nargin
                return
            end
            sL.p = p;
            if nargin > 1
                sL = varassign(sL, varargin);
            end
        end
        %% operational fcns
        function [rad,ang] = LvL(sLA,sLB)
          [rad,ang] = VvV(sLA.vec,sLB.vec);
        end
        
        function d = L_p(sL,ps)
          for i = 1:size(ps,1)
            p = ps(i,:);
            vA = p - sL.p;
            vB = sL.vec;
            d(i) = abs(cross(vA,vB))/norm(vB);
          end
        end
        
        function p = LxL(sLA, sLB) %%SLOW!!!
            p = []; %only consider Point output  
            nA = numel(sLA);
            nB = numel(sLB);
            for i = 1:nA
                for j = 1:nB
                    LA = sLA(i);
                    LB = sLB(j);
                    if isparallel(LA,LB)
                        pnow = LA.p(LB.isin(LA.p),:);
                        if isa(LB,'SegL') 
                            pnow = [pnow; LB.p(LA.isin(LB.p),:); LA.p(LB.isin(LA.p),:)]; 
                        end
                        pnow = uniquei(pnow);
                        p = [p; pnow];
                    else
                        A = [ LA.vec; -LB.vec ];
                        C = LB.p(1,:) - LA.p(1,:);
                        ts = C/A;
                        pA = LA.Lp(ts(1));
                        pB = LB.Lp(ts(2));
                        if eqi(pA,pB) && LA.isin(pA) && LB.isin(pA)
                            p = [p; pA];
                        end
                    end
                end
            end
            p = uniquei(p);
        end
        
        function p = LxE(sL,E)
            p = []; %only consider Point output
            nL = numel(sL);
            nE = size(E,1);
            for i = 1:nL
                for j = 1:nE
                    Equ = E(j,:);
                    if isperpend(sL(i).vec,Equ(1:3)) && isinE(Equ,sL(i).p(1,:)) %L on Equ
                        p = [p; sL(i).p];
                    elseif ~isperpend(sL(i).vec,Equ(1:3))
                        t = (Equ(4) - dot(sL(i).p(1,:),Equ(1:3))) / dot(sL(i).vec, Equ(1:3));
                        pnow = sL(i).Lp(t);
                        if sL(i).isin(pnow)
                            p = [p; pnow];
                        end
                    end
                end
            end
            p = uniquei(p);
        end
        
        function p = Lp(sL,t)
            p = sL.p(1,:)+sL.vec*t;
        end
        
        function h = draw(sL,varargin)
            h = [];
            for i = 1:numel(sL)
                h = [h,varassign(line(sL(i).p(:,1),sL(i).p(:,2),sL(i).p(:,3)),varargin)];
            end
            axis(h(1).Parent,'equal');
        end
        %% is* methods
        function chk = isin(sL,obj)
            if isVec(obj)
                n = size(obj,1);
                chk = false(1,n);
                for i = 1:n
                    p = obj(i,:);
                    vA = sL.p(1,:)-p;
                    vB = sL.p(2,:)-p;
                    chk(i) = dot(vA,vB) <= 0 && isparallel(vA,vB);
                end
            elseif isa(obj, 'SegL')
                n = numel(obj);
                chk = false(1,n);
                for i = 1:n
                    chk(i) = all(isin(obj(i).p));
                end
            end
        end
        
        function chk = isparallel(LA, LB)
            chk = isparallel(LA.vec, LB.vec);
        end
        
        function chk = eq(LA, LB)
            chk = all(ismemberi(LA.p, LB.p));
        end
        
        function L = unique(L)
            i = 1;
            while i <= numel(L) %must check it's number every time
                for j = numel(L) : -1 : i+1
                    if L(i) == L(j)
                        L(j) = [];
                    end
                end
                i = i + 1;
            end
        end
        %% set, get methods
        function sL = set.p(sL,p)
          if isempty(p)
            return
          end
%             assert( size(p,1) == 2, '2 pt generate a sL.' );
            sL.p = p;
            sL.vec = p(2,:) - p(1,:);
            sL.long = norm(sL.vec);
        end
    end
    
    methods (Static)
        function sL = example
            sL = SegL([0 0 0;1 1 1]);
            sL.draw;
        end
        
        function chk = isinchk
            tol = 1.1e-5;
            L = SegL([0 0 0;tol 0 0]);
            p = [0 0 0;...%1
                tol/2 0 0;...%1
                tol 0 0;...%1
                0 tol tol;...
                -tol 0 0;...
                tol/2 tol 0;...
                tol/2 -tol 0;...
                0 tol 0;...
                0 -tol 0];
            chk = L.isin(p);
            L.draw;
            draw(p);
        end
        
        function LxLchk
            tol = Const.Tol;
            La = SegL([0 0 0;1 0 0]);
            Lb = SegL([0.5 0 0;0.5 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0 0 0;0 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0 0 0;-1 0 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([tol 0 0;-1 0 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0.5 -tol 0;0.5 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0.5 tol 0;0.5 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0 tol 0;0 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([0 -tol 0;0 1 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
            Lb = SegL([-tol 0 0;-1 0 0]);
            p = La.LxL(Lb);
            figure;
            La.draw;
            Lb.draw;
            draw(p);
        end
    end
end


