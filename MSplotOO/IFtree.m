classdef IFtree < tree
  properties (Dependent)
    EN
  end
  
  methods
    %% constructors
    function IFT = IFtree(IF,val)
      
      if nargin < 1
        return
      end
      
      if isa(IF, 'tree')
        % Copy constructor
        IFT.Parent = IF.Parent;
        if nargin > 1
          if strcmpi(val, 'clear')
            IFT.Node = cell(numel(IFT.Parent), 1);
          else
            cellval = cell(numel(IFT.Parent), 1);
            for i = 1 : numel(IFT.Parent)
              cellval{i} = val;
            end
            IFT.Node = cellval;
          end
        else
          IFT.Node = IF.Node;
        end
        
      else
        % New object with only root IF
        
        IFT.Node = { IF };
      end
      
    end
    
    function IFT = subIFtree(IFT, MST)
      nID = IFT.BFS;
      for i = 1:IFT.n
        IF = IFT.get(nID(i));
        if ~isempty(IF) && IF.Rank
          Lidx = IFT.subleaves(nID(i));
          Ridx = IFT.subleaves(nID(i+1));
          Lsub = halfsub( IF, MST, Lidx, 'L' );
          Rsub = halfsub( IF, MST, Ridx, 'R' );
          try
            sub = LRinter( IF, Lsub, Rsub );
          catch
            1 
          end
          IF.subIF = sub;
          IFT = IFT.set(nID(i), IF);
        end
      end
    end
    %% operational fcns
    function IFT = chscale(IFT,s)
      for i = 1:IFT.n
        IF = IFT.get(i);
        if isempty(IF)
          continue
        end
        IFN = IF.Equ(1:3);
        IF.p = IF.p*s;
        if ~eqi(IF.Equ(1:3),IFN)
          IF = IF.rev;
        end
        if IF.Rank
          for j = 1:numel(IF.subIF)
            IF.subIF(j).p = IF.subIF(j).p*s;
            if ~eqi(IF.subIF(j).Equ(1:3),IFN)
              IF.subIF(j) = IF.subIF(j).rev;
            end
          end
        end
        IFT = IFT.set(i,IF);
      end
    end
    
    function IFT = Lockchk(IFT)
      mtol = Const.FerroConst.MTol;
      Odr = IFT.BFS;
      for i = 1:IFT.n-1
        IFA = IFT.get(Odr(i));
        IFB = IFT.get(Odr(i+1));
        PA = IFT.getparent(Odr(i));
        PB = IFT.getparent(Odr(i+1));
        if PA~=PB || isempty(IFA) || isempty(IFB)...
            || IFA.Equ(4) < IFB.Equ(4) %is not reverse?
          continue
        end
        R = IFT.subtree(Odr(i)).depth+1;
        mp = mean([IFA.Equ(4); IFB.Equ(4)]);
        IFA = IFA.moveTo(mp-R*mtol);
        IFB = IFB.moveTo(mp+R*mtol);
        IFA.Lock = 'L+';
        IFB.Lock = 'R-';
        IFT = IFT.set(Odr(i),IFA);
        IFT = IFT.set(Odr(i+1),IFB);
      end
    end
    
    function [IFT,chk] = movechk(IFT,Idx,Dis)
      n = length(Idx);
      chk = false(n);
      for i = 1:n
        IF = IFT.get(Idx(i));
        if ~isempty(IF) && Dis(i) &&...
            ~(strcmp(IF.Lock,'L+') && Dis(i) > 0 ...
            || strcmp(IF.Lock,'R-') && Dis(i) < 0)
          chk(i) = true;
          IF.Lock = [];
          IF = IF.move(Dis(i));
          IFT = IFT.set(Idx(i),IF);
        end
      end
    end
    
    function IFT = movetmp(IFT,Idx,Dis)
      for i = 1:length(Idx)
        IF = IFT.get(Idx(i));
        if ~isempty(IF) && Dis(i) &&...
            ~(strcmp(IF.Lock,'L+') && Dis(i) > 0 ...
            || strcmp(IF.Lock,'R-') && Dis(i) < 0)
          IF.Lock = [];
          IF = IF.move(Dis(i));
          IFT = IFT.set(Idx(i),IF);
        end
      end
    end
    
    function nrIF = getnear(IFT, idx)
      nrIF = [];
      IF = IFT.get(idx);
      for i = [1:idx-1 , idx+1:IFT.n] %skip idx
        IFnow = IFT.get(i);
        if ~isempty(IFnow) && IF.isnear(IFnow)
          nrIF = [nrIF, IFnow.subnear(IF)]; %get all subIFnow near IF
        end
      end
    end
    
    function h = draw(IFT, type, varargin)
      if nargin < 2
        type = 'normal';
      end
      if isempty(IFT.get('all'))
        h = [];
      else
        h = IFT.get('all').draw(type, varargin{:});
      end
    end
    
    function EN = getEN_W(IFT,W)
      EN = 0;
      cidx = IFT.getchildren(1);
      for i = 1:numel(cidx)
        IFTtmp = IFT.subtree(cidx(i));
        EN = EN + IFTtmp.EN*W(i);
      end
    end
    %% set, get methods
    function EN = get.EN(IFT)
      EN = 0;
      for i = 1:IFT.n
        IF = IFT.get(i);
        if ~isempty(IF) && isempty(IF.Lock) % assume domain wall merged
          EN = EN + IF.EN;
        end
      end
    end
  end
end

function hsub = halfsub(IF, MST, idx, LR)
  hsub = IFace;
  x = 1;
  for i = 1:length(idx)
    MS = MST.get(idx(i));
    for j = 1:length(MS.BIF)
      BIF = MS.BIF(j);
      if eqi(BIF.Equ, IF.Equ) || eqi(BIF.Equ, IF.rev.Equ)
        hsub(x) = BIF;
        hsub(x).([LR,'Var']) = MS.Var;
        x = x + 1;
        break
      end
    end
  end
end

function sub = LRinter( IF, Lsub, Rsub )
  sub = IFace;
  x = 1;
  for i = 1:length(Lsub)
    for j = 1:length(Rsub)
      temp = Lsub(i).FxF(Rsub(j));
      if ~isa(temp,'Face')
        continue
      end
      if ~eqi(temp.Equ, IF.Equ)
        temp = temp.rev;
      end
      sub(x) = temp;
      sub(x).LVar = Lsub(i).LVar;
      sub(x).RVar = Rsub(j).RVar;
      x = x + 1;
    end
  end
end