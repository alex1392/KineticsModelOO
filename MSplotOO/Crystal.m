classdef Crystal < Part
  properties (Hidden)
    MIFn = [0 1 0];
    fine = 1;
    test = 0;
  end
  
  properties
    Var
    DOF
    CutType
  end
  
  properties (Dependent)
    Rank
  end
  
  methods
    %% constructors
    function obj = Crystal(Vars,DOFs,CutType,varargin)
      if ~nargin
        fprintf( '\nCrystal Constructor:\n' );
        fprintf( 'Use "addMS" to add more MStructs in this crytal, and use "rmMS" to remove.\n' );
        fprintf( 'Use "getTree" to get whole properties of this crystal.\n' );
        return
      else
        obj.Var = Vars;
        obj.DOF = DOFs;
        obj.CutType = CutType;
        obj = varassign(obj,varargin);
      end
    end
    
    function obj = addMS(obj)
      chk = 1;
      while chk == 1
        Var = input('Please input a "Var" of the MS in this crystal......');
        DOF = input('Please input the cooresponding DOF......');
        while log2(length(Var))~=length(DOF)
          DOF = input('The Var and DOF are mismatched.\nPlease input the cooresponding DOF......');
        end
        CutType = input('Please choose the cooresponding CutType......1:COA / 2:EC');
        while CutType ~= 1 && CutType ~= 2
          CutType = input('An unknow input.\nPlease choose the cooresponding CutType......1:COA / 2:EC');
        end
        obj.Var{end+1} = Var;
        obj.DOF{end+1} = DOF;
        obj.CutType{end+1} = CutType;
        chk = input('Continue?  1:yes / otherwise:no......');
      end
    end
    
    function obj = rmMS(obj,i)
      obj.Var(i) = [];
      obj.DOF(i) = [];
    end
    %% operational fcns
    function [obj,MST,IFT] = getTree(obj) %All MSs in line
      n = length(obj.Var);
      if isempty(obj.p)
        obj = defp(obj,n);
      end
      Dp = getDp(obj,n);
      [MST,IFT] = CutMS(obj,Dp,n);
      IFT = IFT.subIFtree(MST);
    end 
    
    function [obj,MST,IFT] = getTree2(obj) %All MSs seperated
      obj = obj.defp;
      MST = MStree(obj);
      IFT = IFtree;
      n = length(obj.Var);
      for i = 1:n
        MS = MStruct(obj.Var{i},obj.DOF{i},obj.CutType{i},'fine',obj.fine,'test',obj.test);
        [MST,IFT] = graftTree(MST,IFT,MS,[]);
      end
    end
    %% set, get methods
    function obj = set.Var(obj,Var)
      if ~iscell(Var)
        Var = {Var};
      end
      obj.Var = Var;
    end
    
    function obj = set.DOF(obj,DOF)
      if ~iscell(DOF)
        DOF = {DOF};
      end
      obj.DOF = DOF;
    end
    
    function obj = set.CutType(obj,CutType)
      if ~iscell(CutType)
        CutType = {CutType};
      end
      obj.CutType = CutType;
    end
    
    function r = get.Rank(obj)
      r = max(cellfun(@length,obj.DOF));
    end
  end
  
  methods (Static)
    function [C, MST, IFT] = example
      p = [1.0000    1.0000    0.7500;
        0.7500    1.0000    0.7500;
        1.0000    0.8750    0.6250;
        1.0000    0.7500    0.6250;
        0.6250    0.8750    0.6250;
        0.1250    0.3750    0.6250;
        0.6250    0.3750    0.6250;
        1.0000    1.0000    0.8750;
        0.7500    1.0000    0.8750;
        1.0000    0.7500    0.8750;
        0.3750    0.6250    0.8750;
        0.8750    0.6250    0.8750];
      C = Crystal({[1 3],[1 2 3 4],[4 3 5 3 1 2 6 2],},...
        {rand,rand(1,2),rand(1,3)},{'COA','COA','COA'},'p',p);
      [C, MST, IFT] = C.getTree;
      MST.draw;
    end
  end
end

function [MST,IFT] = graftTree(MST,IFT,MS,cutIF)
[~,MSTtmp,IFTtmp] = MS.getTree;
IFTtmp = IFTtmp.set(1,cutIF);
MST = MST.graft(1,MSTtmp);
IFT = IFT.graft(1,IFTtmp);
end

function Dp = getDp(obj,n)
d = find(obj.MIFn,1);
pmax = max(obj.p(:,d));
pmin = min(obj.p(:,d));
delta = ( pmax - pmin ) / n;
Dp = pmin +delta : delta : pmax - delta;
end

function [MST, IFT] = CutMS(obj,Dp,n)
MST = MStree(obj);
IFT = IFtree;
Cobj = obj;
for i = 1: n-1
  cfEqu = [ obj.MIFn, Dp(i) ];
  Robj = MStruct(obj.Var{i+1},obj.DOF{i+1},obj.CutType{i+1}, 'fine', obj.fine, 'test', obj.test );
  Lobj = MStruct(obj.Var{i},obj.DOF{i},obj.CutType{i}, 'fine', obj.fine, 'test', obj.test );
  cutIF = IFace([],Lobj.Var,Robj.Var);
  [ Cobj, Lobj, cutIF ] = Cobj.Cut( cfEqu, Robj, Lobj, cutIF );
  [MST, IFT] = graftTree( MST, IFT, Lobj, cutIF );
end
if n == 1
  Cobj = MStruct(obj.Var{1},obj.DOF{1},obj.CutType{1},'fine',obj.fine,'test',obj.test);
end
[MST, IFT] = graftTree(MST,IFT,Cobj,[]);
end