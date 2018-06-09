classdef Face
  properties (Hidden, Constant)
    Tol = Const.FerroConst.Tol;
    MTol = Const.FerroConst.MTol;
  end
  
  properties 
    p
  end
  
  properties %dependent
    Equ
    Edges = SegL;
  end
  
  properties (Hidden) %dependent
    faces
    Area = 0;
    N
  end
  
  properties (Hidden)
    sort = 1;
  end
  
  methods
    %% Constructor
    function obj = Face( p, varargin )
      if ~nargin
        return
      end
      obj.p = p;
      if nargin > 1
        obj = varassign( obj, varargin );
      end
    end
    %% Operational methods
    function [rad,ang] = FvF(FA,FB)
      p = mean(FA.FxF(FB).p);
      pA = mean(FA.p);
      pB = mean(FB.p);
      rad = pvp(p,pA,pB);
      if isnan(rad)
        rad = 0;
      end
      if ~FA.ison(pB)
        rad = 2*pi - rad;
      end
      ang = rad*180/pi;
    end
    
    function p = FxL(F,L)
      p = [];
      nF = numel(F);
      for i = 1:nF
        ptemp = L.LxE(F(i).Equ);
        p = [p; ptemp(F(i).isin(ptemp),:)];
      end
    end
    
    function opt = FxF(FA,FB)
      if isparallel(FA.N,FB.N)
        p = [FA.p( FB.isin(FA.p), :); FB.p( FA.isin(FB.p), :); FA.Edges.LxL(FB.Edges)];
        p = uniquei(p);
        if size(p,1) >= 3
          if isa(FA,'IFace')
            opt = IFace(p,[],[]);
          else
            opt = Face(p);
          end
        elseif size(p,1) == 2
          opt = SegL(p);
        else
          opt = [];
        end
      else
        p = [FA.FxL(FB.Edges); FB.FxL(FA.Edges)];
        p = uniquei(p);
        if size(p,1) == 2
          opt = SegL(p);
        elseif size(p,1) > 2
          save('FxFerr.mat')
          opt = [];
%           error('?')
        else
          opt = [];
        end
      end
    end
    
    function sL = FxE(F,E)
      sL = [];
      n = numel(F);
      for i = 1:n
        Ftmp = F(i);
        L = ExE(Ftmp.Equ, E);
        if ~isa(L,'Line')
          continue
        end
        p = Ftmp.Edges.LxL(L);
        p = uniquei(p);
        if size(p,1) == 2
          sL = [sL, SegL(p)];
        end
      end
    end
    
    function d = F_F(FA,FB)
      d = E_E(FA.Equ,FB.Equ);
    end
    
    function F = rev(F)
      F.Equ = -F.Equ;
    end
    
    function h = draw(obj, type, varargin)
      if ~isa(obj,'IFace')
        obj = IFace(obj.p,[],[]);
      end
      if nargin < 2
        type = 'normal';
      end
      h = [];
      for i = 1:numel(obj)
        F = obj(i);
        if F.Rank
          for j = 1:numel(F.subIF)
            h = [h, F.subIF(j).draw(type, varargin{:})];
          end
        else
          switch type
            case 'normal'
              h = [h, varassign(patch('vertices', F.p,...
                'faces', F.faces, 'facecolor', mean(F.facecolor,1)), varargin)];
            case 'compat'
              h = [h, varassign(patch('vertices', F.p,...
                'faces', F.faces, 'facecolor', repmat(F.compatible,1,3)), varargin)];
            case 'N'
              h = [h, varassign(patch('vertices', F.p,...
                'faces', F.faces, 'facecolor', mean(F.facecolor,1)), varargin)];
              SegL([mean(F.p); mean(F.p) + F.Equ(1:3)]).draw;
            otherwise
              error('?')
          end
        end
      end
      view(h(1).Parent,150,15);
      axis(h(1).Parent,'equal');
    end
    %% is* methods
    function chk = ne(FA, FB)
      chk = ~(FA == FB);
    end
    
    function chk = eq(FA, FB)
      chk = isempty(setdiffi(FA.p, FB.p, 1));
    end
    
    function chk = isnear(FA, FB)
      chk = any(FA.isin(FB.p)) || any(FB.isin(FA.p));
    end
    
    function chk = isin(F, ipt) %%SLOW!!!
      if isVec(ipt)
        chk = logical.empty;
        np = size(ipt,1);
        nfp = size(F.p,1);
        Bp = [F.p; F.p(1,:)];
        for i = 1:np
          p = ipt(i,:);
          if ~isinE(F.Equ, p)
            chk(i) = false;
            continue
          end
          for j = 1:numel(F.Edges)
            if F.Edges(j).isin(p)
              chk(i) = true;
              break
            end
          end
          if length(chk) == i
            continue
          end
          for j = 1:nfp
            vCs(j,:) = cross(Bp(j,:)-p,Bp(j+1,:)-p);
            if any(dot(vCs(1:j-1,:),repmat(vCs(j,:),[j-1 1]),2) < 0)
              chk(i) = false;
              break
            end
          end
          if length(chk) < i
            chk(i) = true;
          end
        end
      elseif isa(ipt, 'SegL') || isa(ipt, 'Face')
        n = numel(ipt);
        chk = false(1,n);
        for i = 1:n
          chk(i) = all(isin(ipt(i).p));
        end
      else
        chk = false;
      end
    end
    
    function chk = ison(F,ipt)
      for i = 1:numel(F)
        chk = isonE(F(i).Equ,ipt);
      end
    end
    %% set, get methods
    function obj = set.p(obj,p) %%SLOW!!!
      assert(isVec(p) && size(p,1) >= 3, 'Invalid input "p".' );
      obj.Equ = p2Equ(p);
      if obj.sort
        ang = p12vp(p,obj.Equ);
        Opt = sortrows([ang,p]);
        obj.p = Opt(:,2:end);
      end
      obj.N = obj.Equ(1:3);
      obj.faces = 1:size(obj.p,1);
      n = size(obj.p,1);
      p = [obj.p; obj.p(1,:)];
      obj.Edges = SegL;
      for i = 1:n
        obj.Edges(i) = SegL(p([i,i+1],:));
      end
      obj.Area = 0;
      for i = 1:n
        obj.Area = obj.Area + cross(p(i,:), p(i+1,:));
      end
      obj.Area = 0.5*norm(obj.Area);
    end
    %% old
    function [rad, ang] = FxFang_old(FA, FB) %!???
      FAN = FA.Equ(1:3);
      FBN = FB.Equ(1:3);
      sL = FA.FxF(FB);
      if isparallel(FAN,FBN)
        if isa(sL, 'SegL')
          rad = pi;
        elseif isa(sL, 'Face')
          rad = 2*pi;
        end
        ang = rad*360/2/pi;
        return
      end
      mp = mean(sL.p);
      CN = cross(FAN, FBN); %common normal
      FAD = cross(CN, FAN)/(norm(CN)); %normalized, lenght = 1
      FBD = cross(CN, FBN)/(norm(CN));
      if ~FA.isin(mp + FAD*Const.Tol)
        FAD = -FAD;
      end %revise
      if ~FB.isin(mp + FBD*Const.Tol)
        FBD = -FBD;
      end
      rad = acos(dot(FAD,FBD));
      bmp = mean(FB.p);
      if ~FA.ison(bmp)
        rad = 2*pi - rad;
      end
      assert(isreal(rad), '?');
      ang = rad*360/2/pi;
    end
  end
  
  methods (Static)
    function obj = example
      obj = Face([1 0 0;0 1 0;0 0 1]);
      obj.draw;
    end
    
    function isinchk
      F = Face([0 0 0;1 0 0;0 1 0;1 1 0;]);
      p = F.p;
      F.isin(p)
      for i = 1:numel(F.Edges)
        F.Edges(i).isin(F.Edges(i).p)
      end
    end
    
    function ang = angchk
      F = Face([0 0 0;0 0 1;1 0 0;1 0 1]);
      FA = Face([0 0 0;0 0 1;1 1 0;1 1 1]);
      FB = Face([0 0 0;0 0 1;0 1 0;0 1 1]);
      FC = Face([0 0 0;0 0 1;-1 1 0;-1 1 1]);
      FD = Face([0 0 0;0 0 1;-1 0 0;-1 0 1]);
      FE = Face([0 0 0;0 0 1;-1 -1 0;-1 -1 1]);
      FF = Face([0 0 0;0 0 1;0 -1 0;0 -1 1]);
      FG = Face([0 0 0;0 0 1;1 -1 0;1 -1 1]);
      F.draw('N')
      FA.draw
      FB.draw
      FC.draw
      FD.draw
      FE.draw
      FF.draw
      FG.draw
      ang = rad2deg([F.FxFang(FA);...
        F.FxFang(FB);...
        F.FxFang(FC);...
        F.FxFang(FD);...
        F.FxFang(FE);...
        F.FxFang(FF);...
        F.FxFang(FG);...
        F.FxFang(F)]);
    end
  end
end
