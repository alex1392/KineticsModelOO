classdef Part
  %% Part  A class represent any region with boundary points and faces.
  
  properties (Hidden,Constant)
    Tol = Const.Tol;
    MTol = Const.MTol;
    vn = Const.vn;
    PlotL = Const.PlotL;
    deffine = Const.deffine;
  end
  
  properties
    p
  end
  
  properties (Hidden)
    %p
    vertices
    faces
    BIF
    Vol
  end
  
  methods
    %% Constructor
    function P = Part(p)
      if ~nargin
        return
      end
      P.p = p;
    end
    
    %% operational fcns
    function P = Merge(Ps)
      % Ps must be a series
      P = Ps(1);
      p = {Ps.p};
      for i = 1:length(p)-1
        ip = intersecti(p{i},p{i+1});
        p{i} = setdiffi(p{i},ip);
        p{i+1} = setdiffi(p{i+1},ip);
      end
      P.p = vertcat(p{:});
    end
    
    function [PR,PL,F] = Cut(P,E,PR,PL,F)
      if nargin < 3
        PR = P;  PL = P;
      end
      if nargin < 5
        F = IFace;
      end
      if isempty(E)
        PR = P;  PL = P;  F = [];
        return
      end
      E = fitE(E,P.p); %10/10 
      F = P.PxE(E,F);
      if ~isa(F,'Face')
        PR = []; PL = []; F = [];
        return
      elseif ~eqi(F.Equ, E)
        F = F.rev;
      end
      p = setdiffi(P.p,F.p); 
      [Rp,Lp] = SeperPoint(p,E);
      if isempty(Rp) || isempty(Lp)
        PR = []; PL = []; F = [];
        return
      end
      PR.p = [Rp;F.p]; 
      PL.p = [Lp;F.p]; 
      
      function E = fitE(E,p)
        for i = 1:size(p,1)
          if isinE(E,p(i,:))
            E(4) = dot(E(1:3),p(i,:)); %move Equ to fit point
            break
          end
        end
      end
      
      function [Rp,Lp] = SeperPoint(p,E)
        n = size(p,1);
        Rp = [];  Lp = [];
        for i = 1:n
          if isonE(E,p(i,:))
            Rp = [Rp; p(i,:)];
          else
            Lp = [Lp; p(i,:)];
          end
        end
      end
    end
    
    function [PR,PL,F] = Cut1010(P,E,PR,PL,F)
      if nargin < 3
        PR = P;  PL = P;
      end
      if nargin < 5
        F = IFace;
      end
      if isempty(E)
        PR = P;  PL = P;  F = [];
        return
      end
      F = P.PxE(E,F);
      if ~isa(F,'Face')
        PR = []; PL = []; F = [];
        return
      elseif ~eqi(F.Equ, E)
        F = F.rev;
      end
      p = setdiffi(P.p,F.p); 
      [Rp,Lp] = SeperPoint(p,E);
      if isempty(Rp) || isempty(Lp)
        PR = []; PL = []; F = [];
        return
      end
      PR.p = [Rp;F.p]; 
      PL.p = [Lp;F.p]; 
      
      function [Rp,Lp] = SeperPoint(p,E)
        n = size(p,1);
        Rp = [];  Lp = [];
        for i = 1:n
          if isonE(E,p(i,:))
            Rp = [Rp; p(i,:)];
          else
            Lp = [Lp; p(i,:)];
          end
        end
      end
    end
    
    function sL = PxL(P, L)
      p = P.BIF.FxL(L);
      p = uniquei(p);
      if isempty(p)
        sL = [];
        return
      elseif size(p,1) == 1
        p(2,:) = L.p(P.isin(L.p));
      end
      assert(size(p,1) == 2,'?');
      sL = SegL(p);
    end
    
    function F = PxE(P,E,F)
      if nargin < 3
        F = Face;
      end
      sL = P.BIF.FxE(E);
      if ~isa(sL,'SegL')
        F = [];
        return
      end
      p = uniquei(vertcat(sL.p));
      if size(p,1) < 3
        F = [];
        return
      end
      F.p = p;
    end
    
    function FB = PxF(P,FA)
      F = P.PxE(FA.Equ);
      FB = FA.FxF(F);
    end
    
    function P = defp(P, n)
      if nargin < 2
        n = 1;
      end
      i = 1;
      p = zeros(8,3);
      for x = [ P.PlotL, -P.PlotL ]/2
        for y = [ P.PlotL, -P.PlotL ]/2*n
          for z = [ P.PlotL, -P.PlotL ]/2
            p(i,:) = [ x, y, z ];
            i = i + 1;
          end
        end
      end
      P.p = p;
    end
    
    function h = draw(P, varargin)
      h = [];
      for i = 1:numel(P)
        temp = patch('vertices', P(i).vertices, 'faces', P(i).faces, 'facecolor', 'w');
        h = [h, varassign(temp, varargin)];
      end
      view(h(1).Parent,150,15);
      axis(h(1).Parent,'equal');
    end
    %% is* methods
    function chk = isin(P, obj)
      if isVec(obj)
        n = size(obj,1);
        chk = false(1,n);
        for i = 1:n
          p = obj(i,:);
          chk(i) = all(P.BIF.ison(p) || P.BIF.isin(p));
        end
      elseif isa(obj, 'SegL') || isa(obj, 'Face')
        n = numel(obj);
        chk = false(1,n);
        for i = 1:n
          chk(i) = P.isin(obj(i).p);
        end
      end
    end
    %% set, get methods
    function P = set.p(P,p)
      assert(isVec(p) && size(p,1) >= 4, 'more than 4 pt generate a Part');
      P.p = p;
      P.vertices = p;
      [faces, fEqu, Vol] = getVolIF(p);
      P.BIF = IFace; %reset
      for f = 1:size(faces,1) % sort faces
        facesN = faces( f, ~isnan(faces(f,:)) );
        facesS = fpSort( p, facesN, fEqu(f,:) );
        faces(f,:) = [ facesS, faces( f, isnan(faces(f,:)) ) ];
        P.BIF(f) = IFace(p(facesS,:),[],[],'sort',0,'Fp',facesS);
      end
      P.faces = faces;
      P.Vol = Vol;
      
      % Set.p Inner Functions
      function [faces, Equ, Vol] = getVolIF(p)
        [faces, Vol] = convhull(p);
        cen = mean(p);
        fn = size(faces,1);
        pn = size(p,1);
        Equ = zeros(fn,4);
        for i = fn : -1 : 1
          kill = 0;
          for j = i+1:size(faces,1)
            if all( ismember( faces(i,1:3), faces(j,:) ))
              faces(i,:) = [];
              Equ(i,:) = [];
              kill = 1;
              break
            end
          end
          if kill
            continue
          end
          Equ(i,:) = p2Equ(p(faces(i,1:3),:));
          if ~isonE(Equ(i,:),cen) % inward normal
            Equ(i,:) = -Equ(i,:);
          end
          k = 4;
          for j = pn : -1 : 1
            if ~ismember(j, faces(i,:)) && isinE(Equ(i,:),p(j,:))
              faces(i,k) = j;
              k = k+1;
            end
          end
        end
        faces(~faces) = nan;
      end
      
      function facesS = fpSort( p, facesN, Equ )
        Ang = p12vp(p(facesN,:),Equ);
        Opt = sortrows([Ang,facesN']);
        facesS = Opt(:,2)';
      end
    end
  end
end