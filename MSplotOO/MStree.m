classdef MStree < tree
  properties (Dependent)
    VarVol
  end
  
  properties (Hidden,Dependent)
    Vol
    EffProp
  end
  
  methods
    %% constructor
    function MST = MStree(MS,val)
      
      if nargin < 1
        return
      end
      
      if isa(MS, 'tree')
        % Copy constructor
        MST.Parent = MS.Parent;
        if nargin > 1
          if strcmpi(val, 'clear')
            MST.Node = cell(numel(MST.Parent), 1);
          else
            cellval = cell(numel(MST.Parent), 1);
            for i = 1 : numel(MST.Parent)
              cellval{i} = val;
            end
            MST.Node = cellval;
          end
        else
          MST.Node = MS.Node;
        end
        
      else
        % New MSTect with only root MS
        MST.Node = {MS};
      end
      
    end
    %% operational fcns
    function MST = chscale(MST,s)
      for i = 1:MST.n
        MS = MST.get(i);
        MS.p = MS.p*s;
        MST = MST.set(i,MS);
      end
    end
    
    function h = draw(MST, varargin)
      idx = MST.findleaves;
      if ~isempty(idx)
        h = MST.get(idx).draw(varargin{:});
      else
        h = MST.get(1).draw(varargin{:});
      end
      view(h(1).Parent,150,15);
      axis(h(1).Parent,'equal');
    end
    
    function Vol = getVol_W(MST,W)
      Vol = sum(MST.getVarVol_W(W));
    end
    
    function VarVol = getVarVol_W(MST,W)
      cidx = MST.getchildren(1);
      VarVol = zeros(1,Const.vn);
      for i = 1:numel(cidx)
        VarVol = VarVol + MST.subtree(cidx(i)).VarVol*W(i);
      end
    end
    
    function EffProp = getEProp_W(MST,W)
      EffProp = getEProp(MST.getVarVol_W(W));
    end
    %% set, get methods
    function Vol = get.Vol(MST)
      Vol = sum(MST.VarVol);
    end
    
    function VarVol = get.VarVol(MST)
      VarVol = zeros(1,Const.FerroConst.vn);
      idxs = MST.findleaves;
      for i = 1:numel(idxs)
        D = MST.get(idxs(i));
        VarVol(D.Var) = VarVol(D.Var) + D.Vol;
      end
    end
    
    function EffProp = get.EffProp(MST)
      EffProp = getEProp(MST.VarVol);
    end
  end
  
end