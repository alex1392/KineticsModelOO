classdef IFace < Face
  properties (Constant,Hidden)
    i_ini = Const.i_ini;
    i_ang = Const.i_ang;
    A = Const.A;
  end

  properties
    LVar
    RVar
    subIF
  end
  
  properties (Hidden) 
    %LRVar
    type
    Rank
    compatible
    facecolor
    %subIF
    EN
    I
    %other
    Lock % L+ / R- / []
    Fp
    dmv
  end
  
  methods
    %% Constructor
    function IF = IFace(p,LVar,RVar,varargin)
      if ~nargin
        return
      end
      if ~isempty(p)
        IF.p = p;
      end
      if ~isempty(LVar)
        IF.LVar = LVar;
      end
      if ~isempty(RVar)
        IF.RVar = RVar;
      end
      if nargin > 3
        IF = varassign(IF,varargin);
      end
    end
    %% Operational methods
    function IF = move(IF,d)
      assert( isscalar(d) && isnumeric(d), 'Invalid input distance.');
      IF.Equ(4) = IF.Equ(4) + d;
    end
    
    function IF = moveTo(IF,d)
      assert( isscalar(d) && isnumeric(d), 'Invalid input distance.');
      IF.Equ(4) = d;
    end
    
    function [sIFA, sIFB] = subnearAB(IFA, IFB)
      sIFA = IFA.subnear(IFB);
      sIFB = IFB.subnear(IFA);
    end
    
    function snIFA = subnear(IFA, IFB)
      %get IFA's subIF which is near IFB
      snIFA = [];
      if ~IFA.isnear(IFB)
        return
      end %IFA must near IFB
      for i = 1:numel(IFA.subIF)
        if IFB.isnear(IFA.subIF(i))
          snIFA = [snIFA, IFA.subIF(i)];
        end
      end
      if isempty(snIFA)
        snIFA = IFA;
      end
    end
    %% set,get methods
    function IF = set.LVar(IF, LVar)
      if isempty(LVar)
        return
      end
%       assert( isrow(LVar) && all(LVar > 0 & LVar <= Part.vn) && ~any(mod(LVar,1))...
%         && isnumeric(LVar) && ~mod(log2(length(LVar)),1),  'Invalid input of "LVar".' );
      IF.LVar = LVar;
    end
    
    function IF = set.RVar(IF, RVar)
      if isempty(RVar)
        return
      end
%       assert( isrow(RVar) && all(RVar > 0 & RVar <= Part.vn) && ~any(mod(RVar,1))...
%         && isnumeric(RVar) && ~mod(log2(length(RVar)),1),  'Invalid input of "RVar".' );
      IF.RVar = RVar;
    end
    
    function r = get.Rank(IF)
      if isempty(IF.LVar) || isempty(IF.RVar)
        r = [];
        return
      end
      L = max([length(IF.RVar),length(IF.LVar)]);
      r = log2(L);
    end
    
    function type = get.type(IF)
      if isempty(IF.Rank)
        type = [];
        return
      end
      if IF.Rank
        type = 'mixed';
      else
        if IF.RVar == IF.LVar
          type = '0';
        elseif ceil(IF.RVar/2) == ceil(IF.LVar/2)
          type = '180';
        else
          type = '90';
        end
      end
    end
    
    function com = get.compatible(IF)
      if isempty(IF.Rank)
        com = [];
        return
      end
      switch IF.type
        case '0'
          com = true;
        case '90'
          N = zeros(1,3);
          N(ceil(IF.LVar/2)) = (-1)^(ceil(mod((IF.LVar+1)/2,1)))/sqrt(2);
          N(ceil(IF.RVar/2)) = (-1)^(ceil(mod((IF.RVar+1)/2,1)))/sqrt(2);
          if eqi(N, IF.Equ(1:3)) || eqi(N, -IF.Equ(1:3))
            com = true;
          else
            com= false;
          end
        case '180'
          if (IF.LVar==1 || IF.LVar==2) && eqi(IF.Equ(1), 0) ||...
              (IF.LVar==3 || IF.LVar==4) && eqi(IF.Equ(2), 0) ||...
              (IF.LVar==5 || IF.LVar==6) && eqi(IF.Equ(3), 0)
            com = true;
          else
            com = false;
          end
        case 'mixed'
          com = 'COA';
      end
      
    end
    
    function fc = get.facecolor(IF)
      fc = ones(2,3);
      if IF.Rank
        return
      end
      if ~isempty(IF.LVar)
        fc(1,:) = Const.varcolor(IF.LVar,:);
      end
      if ~isempty(IF.RVar)
        fc(2,:) = Const.varcolor(IF.RVar,:);
      end
    end
    
    function IF = set.subIF(IF,subIF)
      %assert(eqi(sum([subIF.Area]), sum([IF.Area])),'!');
      IF.subIF = subIF;
    end
    
    function EN = get.EN(IF)
      EN = [];
      if isempty(IF.type)
        return
      end
      switch IF.type
        case {'0','90','180'}
          EN = getEN(IF);
%           EN = IF.EN;
        case 'mixed'
          EN = sum(getEN(IF.subIF));
%           EN = sum(IF.subIF.EN)
      end
    end
    
    function I = get.I(IF)
      I = [];
      if isempty(IF.type)
        return
      end
      switch IF.type
        case {'0','90','180'}
          I = getI(IF);
        case 'mixed'
          I = sum(getI(IF.subIF) .* [IF.subIF.Area] ./ IF.Area);
      end
    end
    
  end
end


function EN = getEN(IF)
% matlab can't call get method recursively...
EN = [];
for i = 1:numel(IF)
  deltaP = Const.MProp(IF(i).LVar).PR - Const.MProp(IF(i).RVar).PR;
  EN = [EN, norm(deltaP)^2*IF(i).A(IF(i).compatible+1)*IF(i).Area/2];
end
end

function I = getI(IF)
% matlab can't call get method recursively...
I = [];
for i = 1:numel(IF)
  switch IF(i).type
    case '0'
      I = [I, Const.i_ini*Const.i_ang(1)];
    case '90'
      I = [I, Const.i_ini*Const.i_ang(2)];
    case '180'
      I = [I, Const.i_ini*Const.i_ang(3)];
  end
end
end