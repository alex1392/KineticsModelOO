classdef MStruct < Part
  properties (Hidden)
    fine = 1;
    test = 0;
  end
  
  properties
    Var
    DOF
    CutType = 'COA';
  end
  
  properties (Dependent)
    Rank
    facecolor
  end
  
  methods
    %% constructor
    function MS = MStruct(Var,DOF,CutType,varargin)
      if ~nargin
        fprintf('\nMStruct Constructor:\n');
        fprintf('Please assign "Var" and "DOF" to this MS, then use "getTree" method to get its whole properties.\n');
      end
      if nargin
        MS.Var = Var;
      end
      if nargin >= 2
        MS.DOF = DOF;
      end
      if nargin >= 3
        MS.CutType = CutType;
      end
      if nargin >= 4
        MS = varassign( MS, varargin );
      end
    end
    %% operational fcns
    function [MS,MST,IFT] = getTree(MS)
      if isempty(MS.p)
        MS = defp(MS);
      end
      if ~MS.Rank
        MST = MStree(MS);
        IFT = IFtree(MST,[]);
        return
      end
      assert(MS.Rank == length(MS.DOF), 'Input "Var" and "DOF" are mismatched');
      [Stree,Badness,~] = BuildTree(MS);
      assert(Badness < Const.FerroConst.Tol, 'This MS is not compatible.');
      switch MS.CutType
        case 'COA'
          MST = MakeEmptyTree(MS.Rank, MStree(MS), 1);
          IFT = IFtree(MST, []);
          PFIdx = 2;
          COACut(MS,Stree);
        case 'EC'
          MST = MStree(MS);
          IFT = IFtree();
          
        case '' %Define your cut type here
          
      end
      IFT = IFT.subIFtree(MST);
      
      % getTree Inner Functions
      function ECCut(MST,Stree,Rank)
        if ~Rank
          return
        end
        StartP = mean(MS.p);
        for Node = MST.Node
          if Node.Rank == Rank && Node.isin(StartP)
            break;
          end
        end
        
        
        
      end
      
      function MST = MakeEmptyTree(Rank,MST,Idx)
        if ~Rank
          return
        end
        if strcmp(MS.fine,'rand')
          fine = randi(4);
        elseif strcmp(MS.fine,'input')
          fine = input('input fine\n');
        elseif length(MS.fine) >= Rank
          fine = MS.fine(Rank);
        else
          fine = MS.fine(1);
        end
        n = 2*fine+1;
        x = zeros(1,n);
        for i = 1:n
          [MST,x(i)] = MST.addnode(Idx,[]);
        end
        for i = 1:n
          MST = MakeEmptyTree(Rank-1,MST,x(i));
        end
      end
            
      function COACut(MS,Stree)
        if ~MS.Rank
          return
        end
        [LST,RST] = SetTree(MS,Stree);
        CpartR = CutRankCOA(MS,Stree);
        for CID = 1:length(CpartR)
          if mod(CID,2) %odd
            COACut(CpartR(CID),LST);
          else %even
            COACut(CpartR(CID),RST);
          end
        end
      end
      
      function [LST,RST] = SetTree(part,Stree)
        RST = [];  LST = [];
        for r = 1:part.Rank
          RST = [RST,Stree(2^r : 2^(r+1) - 2^(r-1) - 1)];
          LST = [LST,Stree(2^(r+1) - 2^(r-1) : 2^(r+1) - 1)];
        end
      end
      
      function CpartR = CutRankCOA(parent,Stree)
        IFN = Stree(1).IfNormal(:,1);
        Dp = findDpCOA(parent,IFN);
        CpartR = paraCutCOA(parent,IFN,Dp);
      end
      
      function Dp = findDpCOA(parent,cfNormal)
        p = parent.p;
        n = size(p,1);
        d = zeros(n,1);
        for i = 1:n
          d(i) = dot(p(i,:), cfNormal); % faceEqu constant
        end
%         if strcmp(parent.fine,'rand')
%           fine = randi(3);
%         elseif length(parent.fine) >= parent.Rank
%           fine = parent.fine(parent.Rank);
%         else
%           fine = parent.fine(1);
%         end
fine = (numel(MST.getchildren(MST.Parent(PFIdx)))-1)/2;
        deltaD = (max(d) - min(d)) / fine;
        Dp = min(d) : deltaD : max(d);
        deltaD2 = 0.5*parent.DOF(end)*deltaD;
        for i = 2:2: 2*(fine - 1)
          Dp = [ Dp(1:i-1), Dp(i)-deltaD2, Dp(i)+deltaD2, Dp(i+1:end) ];
        end
        Dp = [ Dp(1) + deltaD2, Dp(2:end-1), Dp(end) - deltaD2 ];
      end
      
      function CpartR = paraCutCOA(parent,IFN,Dp)
        CpartR = [];
        Cobj = parent;
        for i = 1: length(Dp)
          cfEqu = [IFN', Dp(i)];
          [Robj,Lobj,CIF] = PreSet(parent,i);
          [Robj,Lobj,CIF] = Cobj.Cut(cfEqu,Robj,Lobj,CIF);
          CpartR = [CpartR,Lobj];
          MST = MST.set(PFIdx,Lobj); %Lobj.draw
          IFT = IFT.set(PFIdx,CIF);
          PFIdx = PFIdx + 1;
          Cobj = Robj;
        end
        CpartR = [CpartR,Robj];
        MST = MST.set( PFIdx, Robj ); %Robj.draw
        PFIdx = PFIdx + 1;
      end
      
      function [Robj,Lobj,CIF] = PreSet(parent,cutID)
        Robj = parent;
        Lobj = Robj;
        Robj.DOF = parent.DOF(1:end-1);
        Lobj.DOF = parent.DOF(1:end-1);
        RVar = parent.Var( (length( parent.Var )/2)+1 : end );
        LVar = parent.Var( 1 : (length( parent.Var )/2) );
        if mod(cutID,2)
          Robj.Var = RVar;
          Lobj.Var = LVar;
        else
          Robj.Var = LVar;
          Lobj.Var = RVar;
        end
        CIF = IFace([],Lobj.Var,Robj.Var);
      end
    end
    
    function h = draw(MS, varargin)
      h = [];
      for i = 1:numel(MS)
        temp = patch('vertices', MS(i).vertices, 'faces', MS(i).faces, 'facecolor', MS(i).facecolor);
        h = [h, varassign(temp, varargin)];
      end
      view(h(1).Parent,150,15);
      axis(h(1).Parent,'equal');
      %             drawnow
    end
    %% set, get methods
    function fc = get.facecolor(MS)
      if MS.Rank
        fc = 'w';
      else
        fc = Const.FerroConst.varcolor(MS.Var,:);
      end
    end
    
    function r = get.Rank(MS)
      r = length(MS.DOF);
    end
  end
  
  methods (Static)
    function [MS,MST,IFT] = example
      %             p = [1.0000    1.0000    0.7500;
      %                 0.7500    1.0000    0.7500;
      %                 1.0000    0.8750    0.6250;
      %                 1.0000    0.7500    0.6250;
      %                 0.6250    0.8750    0.6250;
      %                 0.1250    0.3750    0.6250;
      %                 0.6250    0.3750    0.6250;
      %                 1.0000    1.0000    0.8750;
      %                 0.7500    1.0000    0.8750;
      %                 1.0000    0.7500    0.8750;
      %                 0.3750    0.6250    0.8750;
      %                 0.8750    0.6250    0.8750];
      MS = MStruct([4 3 5 3 1 2 6 2],[0.5 0.5 0.5],'COA','fine',1);
      [MS, MST, IFT] = MS.getTree;
      figure;
      MST.draw;
      figure;
      IFT.draw;
    end
  end
end