function [h,F] = draw(obj,sz,varargin)
%%
hold(gca,'on');
view(150,15);
axis('equal');
F = [];
if isVec(obj)
  if nargin < 2
    sz = 36;
  end
  h = scatter3(obj(:,1),obj(:,2),obj(:,3),sz,'filled',varargin{:});
elseif isEqu(obj)
  if nargin < 2
    sz = 2;
  end
  [h,F] = drawEqu(obj,sz);
end

%%
  function [h,F] = drawEqu(obj,Size)
    %%
    h = [];
    for i = 1:size(obj,1)
      Equ = obj(i,:);
      vA = Equ(1:3);
      pA = ( Equ(4) / sum(vA.^2) )*vA;
      
      %vB =  (vA \ Equ(4))' - pA;
      vB = [0 0 0];
      while ~any(vB)
        pB = randPtOnEqu(Equ);
        vB = pB - pA;
      end
      
      vC = cross(vA,vB);
      uB = vB/norm(vB)*Size;
      uC = vC/norm(vC)*Size;
      F = IFace([ pA+uB; pA+uC; pA-uB; pA - uC ],[],[]);
      h = [h,F.draw('normal',varargin{:})];
    end
  end

  function p = randPtOnEqu(E)
    assert(any(E),'randPtOnEqu error!');
    if E(3) ~= 0
      x = rand;
      y = rand;
      z = E(4)-(E(1)*x + E(2)*y)/E(3);
    elseif E(2) ~= 0
      x = rand;
      z = rand;
      y = E(4)-(E(1)*x + E(3)*z)/E(2);
    elseif E(1) ~= 0
      y = rand;
      z = rand;
      x = E(4)-(E(3)*z + E(2)*y)/E(1);
    end
    p = [x,y,z];
  end
end