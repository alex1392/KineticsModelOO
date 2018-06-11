function h = draw(obj,sz,varargin)
    %%
    hold(gca,'on');
    view(150,15);
    axis('equal');
    if isVec(obj)
      if nargin < 2
        sz = 36;
      end
        h = scatter3(obj(:,1),obj(:,2),obj(:,3),sz,'filled',varargin{:});
    elseif isEqu(obj)
      if nargin < 2
        sz = 2;
      end
        h = drawEqu(obj,sz);
    end

    %%
    function h = drawEqu(obj,Size)
        %%
        h = [];
        for i = 1:size(obj,1)
            Equ = obj(i,:);
            vA = Equ(1:3);
            pA = ( Equ(4) / sum(vA.^2) )*vA;
            vB =  (vA \ Equ(4))' - pA;
            vC = cross(vA,vB);
            uB = vB/norm(vB)*Size;
            uC = vC/norm(vC)*Size;
            F = IFace([ pA+uB; pA+uC; pA-uB; pA - uC ],[],[]);
            h = [h,F.draw('normal',varargin{:})];
        end
    end

end