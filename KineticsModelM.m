function KineticsModel()
%% KINETICSMODELOO  Simulate the hysteresis behavior of ferroelectric crystals, object-oriented version.
%% Initialize Matlab Environment
clear global; close all;
addpath(genpath(pwd))
%% Initial State
varR0 = []; varR1 = [];
dofR0 = []; dofR1 = [];
load('vardof')
if ~nargin
  var = {[1 2],[1 3],[2 3]};
  dof = {0.5,0.5,0.5};
  type = repmat({'COA'},[1 length(var)]);
end
obj = Crystal(var,dof,type);
CN = numel(obj.Var);
%% Load System
LoadDeg = 1e6; % Mega
EVec = [1;0;0];  % MV
sigma = [0;0;0;0;0;0];  % MPa
DirecE = find(EVec,1);
if isempty(DirecE)
  DirecE = 3;
end
dir = num2str(DirecE);
%% Options
IFMethod = 'N'; % 'N', 'S'
MSWeight = 1;
opt = '-p \.jpg'; % '-d' = data('\*'), '-p' = plot('\.fig','\.jpg'), '-m' = movie('\.mp4','\.gif')
IFtype = 'compat'; %'normal', 'compat'
IFMerge = 0;
%% Parameters
sEF = 1e6;
fps = 1; %frames per secend
ScaleL=1e-2; %true scale (m)
Lfrac = ScaleL / obj.PlotL;  %scale fraction
dt = 1e-3; % spacing time
ti = 0; % starting time
tf = 1.25; % ending time
T = 0; %true time
%% Figure Setting
ax = [];
hT = [];
hRT = [];
fcharts = [];
fSC = [];
fIF = [];
%% Data Output
if CN == 1
  VarStr = num2str( obj.Var{1} );
elseif CN > 1
  VarStr = [ num2str(CN),'SC' ];
end
if isempty(strfind(pwd,VarStr))
  mkdir(VarStr); cd(VarStr)
end
fName = [VarStr,'_',num2str(sigma(DirecE)),'MPa' ];
if strfind( opt, '-d' )
  fid = fopen([ fName, '.xls' ], 'wt' );
  fprintf( fid, 'KineticsModel\n' );
  
  fprintf( fid, '\n%%Initial State\n' );
  fprintf( fid, 'Var:\t' );
  for CNID = 1:CN
    fprintf( fid, '%s\t', num2str( obj.Var{CNID} ) );
  end
  fprintf( fid, '\nDOF:\t' );
  for CNID = 1:CN
    fprintf( fid, '%s\t', num2str( obj.DOF{CNID} ) );
  end
  
  fprintf( fid, '\n%%Load System\n' );
  fprintf( fid, 'Load Degree:\t%d\n', LoadDeg );
  fprintf( fid, 'Electric field:\t%s\n', num2str( EVec' ) );
  fprintf( fid, 'Stress tensor:\t%s\n', num2str( sigma' ) );
  
  fprintf( fid, '\n%%Options\n' );
  fprintf( fid, 'IFMethod:\t%s\n', IFMethod );
  
  fprintf( fid, '\n%%Parameters\n' );
  fprintf( fid, 'dt:\t%d\n', dt );
  fprintf( fid, 'Ps:\t%d\n', obj.EffProp.Ps );
  fprintf( fid, 'epss:\t%d\n', obj.EffProp.epss );
  fprintf( fid, 'fine:\t%d\n', obj.fine );
  fprintf( fid, 'deffine:\t%s\n', num2str( obj.deffine ) );
  fprintf( fid, 'inertia:\t%d\n', Face.i_ini );
  fprintf( fid, 'i_ang:\t%s\n', num2str( Face.i_ang ) );
  fprintf( fid, 'fps:\t%d\n', fps );
  
  fprintf( fid, '\n%%Simulation Datas\n' );
  fprintf( fid, '\nt\tEF\tD\teps\tVar1\tVar2\tVar3\tVar4\tVar5\tVar6\t\n' );
end
if strfind( opt, '-m' )
  frames = struct( 'cdata', [], 'colormap', [] );
  if strfind( opt, '\.mp4' )
    movie = VideoWriter(fName,'MPEG-4');
    movie.Quality = 100;
    movie.FrameRate = 15;
  elseif strfind( opt, '\.gif' )
  end
end
%% Simulation Starts
if ~MSWeight
  [obj, MST, IFT] = obj.getTree;
else
  [obj, MST, IFT] = obj.getTree2;
end
axLim = [min(obj.p(:,1)),max(obj.p(:,1)),...
      min(obj.p(:,2)), max(obj.p(:,2)),...
      min(obj.p(:,3)), max(obj.p(:,3))];
axLim = axLim*CN^(1/3) + [-1 1 -1 1 -1 1]*0.1;
picN = 1;
t = ti;
while t < tf
  stress = sigma*LoadDeg; %load mode
  EF = EVec*sin(2*pi*t)*LoadDeg;
  if t == ti
    data = setfig(MST,IFT,EF,stress);
    V = MST.Vol;
  elseif abs(EF(DirecE)) < sEF
    tic
    IFM = judgeIF(MST,IFT,EF,stress);
    [MST,IFT] = IFmoveAll(MST,IFT,IFM);
    dW = judgeWeight(MST,IFT,EF,stress);
    [MST,IFT] = moveWeight(MST,IFT,dW);
    if ~eqi(MST.Vol,V)
      1
    end
    data = replot(data,MST,IFT,EF,stress);
  else
    t = t + dt;
    continue
  end
  if eqi(MST.EffProp.epsR(DirecE),MST.EffProp.epss) && sEF == LoadDeg
      sEF = 2*EF(DirecE);
  end
  if ~isempty(opt)
    if strfind( opt, '-d' )
      if ~( strfind(opt, '\*') && t < 0.25 )
        VarVol = arrayfun(@(x) {x}, MST.VarVol);
        fprintf( fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',t ,...
          EF(DirecE)/LoadDeg*10,MST.EffProp.PR(DirecE) ,MST.EffProp.epsR(DirecE)*100,VarVol{:});
      end
    end
    if strfind( opt, '-p' )
      if strfind(opt, '\.fig')
        saveas([fSC, fIF, fcharts], [fName, '-', num2str(picN), '.fig']);
      elseif strfind(opt, '\.jpg')
        saveas(fIF, [fName, '-IF', num2str(picN), '.jpg']);
        saveas(fSC, [fName, '-SC', num2str(picN), '.jpg']);
        saveas(fcharts, [fName, '-charts', num2str(picN), '.jpg']);
      end
    end
    if strfind( opt, '-m' )
      frames(picN) = getframe(fSC);
    end
  end
  pause(dt);
  picN = picN + 1;
  t = t + dt;
  if ~mod(picN,10)
    save([fName,'_',num2str(picN)])
  end
end
if strfind( opt, '-d' )
  fclose(fid);
end
if strfind( opt, '-m' )
  if strfind( opt, '\.mp4' )
    open(movie);
    writeVideo(movie, frames);
    close(movie);
  elseif strfind( opt, '\.gif' )
    for fid = 1:numel(frames)
      [image,~]=frame2im(frames(fid));
      [im,map]=rgb2ind(image,128);
      if fid==1
        imwrite(im,map,[ fName, '.gif' ],'gif','writeMode','overwrite','delaytime',0.1,'loopcount',inf);
      else
        imwrite(im,map,[ fName, '.gif' ],'gif','writeMode','append','delaytime',0.1);
      end
    end
  end
end
%% KineticsModelOO Inner Functions
  function IFM = judgeIF(MST,IFT,EF,stress)
    n = IFT.nnodes;
    d = Const.MTol/2;
    IFM = zeros(1,n);
    if ~MSWeight
      switch IFMethod
        case 'N' %energy + draw
          for i = 1: n
            if isempty(IFT.get(i)) || IFT.depthtree.get(i) == 1 %SC domain can't move
              continue
            end
            IF = IFT.get(i);
            [MStempR,IFtempR] = IFmove(MST,IFT,i,d);
            [MStempL,IFtempL] = IFmove(MST,IFT,i,-d);
            GR = getG(MStempR.VarVol,EF,stress) + IFtempR.EN;
            GL = getG(MStempL.VarVol,EF,stress) + IFtempL.EN;
            g = (GR - GL) * MST.Vol * Lfrac^3 / (2*d);
            IFM(i) = -g/(IF.I*IF.Area*Lfrac^2)*dt*Lfrac;
            if abs(IFM(i)) > 0.05
              %             save(['IFMerr-',num2str(t),'.mat']);
              IFM(i) = IFM(i)/abs(IFM(i))*0.05;
            end
          end
        case 'S' %Vfrac simplify
          for i = 1:n
            if isempty(IFT.get(i)) || IFT.depthtree.get(i) == 1 %SC domain can't move
              continue
            end
            IF = IFT.get(i);
            bIFs = IFT.getnear(i);
            [RVarVol, LVarVol] = moveVF(MST,IF,d);
            [RwallEN, LwallEN] = movewall(IFT, IF, bIFs, d);
            GR = getG(RVarVol,EF,stress) + RwallEN;
            GL = getG(LVarVol,EF,stress) + LwallEN;
            g = (GR - GL)*MST.Vol*Lfrac^3 / (2*d);
            IFM(i) = -g/(IF.I*IF.Area*Lfrac^2)*dt*Lfrac;
            if abs(IFM(i)) > 0.05
              %             save(['IFMerr-',num2str(t),'.mat']);
              IFM(i) = IFM(i)/abs(IFM(i))*0.05;
            end
          end
      end
    else
      switch IFMethod
        case 'N' %energy + draw
          for i = 1: n
            if isempty(IFT.get(i)) || IFT.depthtree.get(i) == 1 %SC domain can't move
              continue
            end
            IF = IFT.get(i);
            [MStempR,IFtempR] = IFmove(MST,IFT,i,d);
            [MStempL,IFtempL] = IFmove(MST,IFT,i,-d);
            GR = getG(MStempR.VarVol,EF,stress) + IFtempR.EN;
            GL = getG(MStempL.VarVol,EF,stress) + IFtempL.EN;
            g = (GR - GL) * MST.Vol * Lfrac^3 / (2*d);
            IFM(i) = -g/(IF.I*IF.Area*Lfrac^2)*dt*Lfrac;
            if abs(IFM(i)) > 0.05
              %             save(['IFMerr-',num2str(t),'.mat']);
              IFM(i) = IFM(i)/abs(IFM(i))*0.05;
            end
          end
        case 'S' 
          for i = 1:n
            if isempty(IFT.get(i)) || IFT.depthtree.get(i) == 1 %SC domain can't move
              continue
            end
            IF = IFT.get(i);
            cID = IFT.getancestors(i); 
            cID = cID(end-1); %end = SC, end-1 = MS
            cIFT = IFT.subtree(cID);
            bIFs = cIFT.getnear(i-cID+1);
            [RVarVol, LVarVol] = moveVF(MST,IF,d);
            [RwallEN, LwallEN] = movewall(IFT, IF, bIFs, d);
            GR = getG(RVarVol,EF,stress) + RwallEN;
            GL = getG(LVarVol,EF,stress) + LwallEN;
            g = (GR - GL)*MST.Vol*Lfrac^3 / (2*d);
            IFM(i) = -g/(IF.I*IF.Area*Lfrac^2)*dt*Lfrac;
            if abs(IFM(i)) > 0.05
              %             save(['IFMerr-',num2str(t),'.mat']);
              IFM(i) = IFM(i)/abs(IFM(i))*0.05;
            end
          end
      end
    end
  end

  function [MST,IFT] = IFmoveAll(MST,IFT,d)
    Idx = 1:IFT.nnodes;
    IFT = IFT.movechk(Idx,d);
    if ~IFMerge
      IFT = IFT.Lockchk;
    else
      [MST,IFT] = Mergechk(MST,IFT);
    end
    [MST,IFT] = ReCut(MST,IFT);
    if IFMerge
      [MST,IFT] = Parentchk(MST,IFT);
      [MST,IFT] = LRVarchk(MST,IFT);
    end
    IFT = IFT.subIFtree(MST);
  end
  
  function dW = judgeWeight(MST,IFT,EF,stress)
    d = Const.MTol/2;
    cidx = MST.getchildren(1);
    Weight = [MST.get(cidx).Vol];
    wallI = Const.i_ini;
    wallA = (MST.Vol/CN)^(2/3);
    dV = d*wallA;
    VF = MST.VarVol;
    for i = 1:CN
      for j = i+1:CN
        MSi = MST.subtree(cidx(i));
        MSj = MST.subtree(cidx(j));
        VFi = MSi.VarVol;
        VFj = MSj.VarVol;
        dVFi = dV*VFi/sum(VFi);
        dVFj = dV*VFj/sum(VFj);
        RVFi = VFi + dVFi; RVFj = VFj - dVFj;
        LVFi = VFi - dVFi; LVFj = VFj + dVFj;
        RVFi(RVFi<0) = 0; RVFj(RVFj<0) = 0;
        LVFi(LVFi<0) = 0; LVFj(LVFj<0) = 0;
        RVarVol = (RVFi-VFi) + (RVFj-VFj) + VF;
        LVarVol = (LVFi-VFi) + (LVFj-VFj) + VF;
        ENi = IFT.subtree(cidx(i)).EN;
        ENj = IFT.subtree(cidx(j)).EN;
        RwallEN = (RVFi/VFi)^(2/3)*ENi + (RVFj/VFj)^(2/3)*ENj;
        LwallEN = (LVFi/VFi)^(2/3)*ENi + (LVFj/VFj)^(2/3)*ENj;
        GR = getG(RVarVol,EF,stress) + RwallEN;
        GL = getG(LVarVol,EF,stress) + LwallEN;
        g = (GR - GL)*MST.Vol*Lfrac^3 / (2*d);
        moved = -g/(wallI*wallA*Lfrac^2)*dt*Lfrac;
        deltaV(i,j) = moved*wallA;
        deltaV(j,i) = -moved*wallA;
      end
    end
    out = 0;
    while out < CN
      out = 0;
      for i = 1:CN
        tmp = Weight(i) + sum(deltaV(i,:));
        if tmp < 0
          eater = find(deltaV(i,:) < 0);
          deltaV(eater,i) = deltaV(eater,i) + (tmp-Const.Tol)/numel(eater); %Wi < 0
          deltaV(i,eater) = deltaV(i,eater) - (tmp-Const.Tol)/numel(eater); %Wi < 0
        else
          out = out + 1;
        end
      end
    end
    Weight = Weight + sum(deltaV,2)';
    if any(Weight < 0)
      1
    end
    dW = Weight - [MST.get(cidx).Vol];
  end
  
  function [MST,IFT] = moveWeight(MST,IFT,dW)
    cidx = MST.getchildren(1);
    for i = 1:CN 
      sMST = MST.subtree(cidx(i));
      factor = ((sMST.Vol+dW(i))/sMST.Vol)^(1/3);
      try
      chMST = sMST.chscale(factor);
      catch
        1
      end
      MST = MST.setsub(cidx(i),chMST);
      chIFT = IFT.subtree(cidx(i)).chscale(factor);
      IFT = IFT.setsub(cidx(i),chIFT);
    end
  end

  function [MST,IFT] = Parentchk(MST,IFT)
    Odr = MST.BFS;
    for i = MST.n:-1:2
      MS = MST.get(Odr(i));
      cidx = MST.getchildren(Odr(i));
      if isempty(cidx)
        continue
      elseif length(cidx) == 1
        cMS = MST.get(cidx);
        MS.Var = cMS.Var;
        MS.DOF = cMS.DOF;
        MST = MST.set(Odr(i),MS);
        MST = MST.chop(cidx);
        IFT = IFT.chop(cidx);
      end
    end
  end

  function [MST,IFT] = LRVarchk(MST,IFT)
    Odr = IFT.BFS;
    for i = IFT.n:-1:1
      IF = IFT.get(Odr(i));
      if isempty(IF)
        continue
      end
      IF.LVar = MST.get(Odr(i)).Var;
      IF.RVar = MST.get(Odr(i+1)).Var;
    end
  end

  function [MST,IFT] = IFmove(MST,IFT,Idx,d)
    IFT = IFT.movechk(Idx,d);
    if ~IFMerge
      IFT = IFT.Lockchk;
    else
      [MST,IFT] = Mergechk(MST,IFT);
    end
    [MST, IFT] = ReCutS(MST,IFT,Idx);
    if IFMerge
      [MST,IFT] = Parentchk(MST,IFT);
      [MST,IFT] = LRVarchk(MST,IFT);
    end
    IFT = IFT.subIFtree(MST);
  end

  function [MST,IFT] = Mergechk(MST,IFT)
    Odr = IFT.BFS;
    i = IFT.n - 2;
    while i > 0
      if IFT.getparent(Odr(i)) == IFT.getparent(Odr(i+1)) &&...
          ~isempty(IFT.get(Odr(i))) &&...
          ~isempty(IFT.get(Odr(i+1))) &&...
          IFT.get(Odr(i)).Equ(4) > IFT.get(Odr(i+1)).Equ(4) %is reversed?
        MSs = MST.get(Odr(i:i+2));
        MS = MSs.Merge;
        MST = MST.set(Odr(i+2),MS);
        IFT = IFT.chop(Odr(i+1)); %backward elimination
        IFT = IFT.chop(Odr(i));
        MST = MST.chop(Odr(i+1));
        MST = MST.chop(Odr(i));
        i = i - 1;
      end
      i = i - 1;
    end
  end

  function [MST,IFT] = ReCut(MST,IFT)
    idx = MST.BFS;
    if MSWeight
      idx(1) = [];
    end
    for i = idx
      if i > IFT.n 
        return
      end
      [MST, IFT] = parallelCut(MST, IFT, i);
    end
  end

  function [MST,IFT] = parallelCut(MST,IFT,Idx)
    Cobj = MST.get(Idx);
    cIdx = MST.getchildren(Idx);
    n = numel(cIdx);
    kill = 0;
    for i = n-1:-1:1
      IF = IFT.get(cIdx(i));
      Rin = MST.get(cIdx(i+1));
      Lin = MST.get(cIdx(i));
      [Rout,Lout,CutIF] = Cobj.Cut(IF.Equ,Rin,Lin,IF);
      if isempty(Rout) %out
        if ~IFMerge
          Rank = IFT.subtree(i).depth+1;
          [Rout,Lout,CutIF] = IFBR(Cobj,IF,Rin,Lin,Rank,1);
          MST = MST.set( cIdx(i+1), Rout );
          IFT = IFT.set( cIdx(i), CutIF );
          Cobj = Lout;
        else
          if i == n-1 %R
            MSs = MST.get(cIdx(i:i+1));
            MS = MSs.Merge;
            MST = MST.set(cIdx(i),MS);
            IFT = IFT.chop(cIdx(i));
            MST = MST.chop(cIdx(i+1));
          elseif i == 1 %L
            kill = 1;
            MSs = MST.get(cIdx(i:i+1));
            MS = MSs.Merge;
            MST = MST.set(cIdx(i+1),MS);
            IFT = IFT.chop(cIdx(i));
            MST = MST.chop(cIdx(i));
          end
        end
      else
        MST = MST.set( cIdx(i+1), Rout );
        IFT = IFT.set( cIdx(i), CutIF );
        Cobj = Lout;
      end
    end
    if n > 1 && ~kill
      MST = MST.set( cIdx(i), Lout );
    end
  end

  function [Rpart,Lpart,CutIF] = IFBR(cpart,IF,Rpart,Lpart,Rank,fail)
    % IFBR  Move the interface back to the edge of the region.
    cfEqu = IF.Equ;
    n = length( cpart.p );
    Dis = zeros( n, 1 );
    for i = 1:n
      Dis(i) = dot( cfEqu(1:3), cpart.p(i,:) ) - cfEqu(4);
    end
    [m,j] = min(abs(Dis));
    [M,k] = max(abs(Dis));
    if M - m > 10*cpart.MTol
      delta =  cpart.MTol*Rank; %delta = (M-m)*0.01;
    else
      delta = (M-m)*0.5;
    end
    Dm = dot( cfEqu(1:3), cpart.p(j,:) );
    DM = dot( cfEqu(1:3), cpart.p(k,:) );
    if Dm < DM
      Dm = Dm + delta;
    elseif Dm > DM
      Dm = Dm - delta;
    end
    cfNew = [ cfEqu(1:3), Dm ];
    [ Rtemp, Ltemp, CutIF ] = cpart.Cut( cfNew, Rpart, Lpart,IF);
    if fail < 5
      if isempty(Rtemp)
        fail = fail + 1;
        [ Rpart, Lpart, CutIF ] = IFBR( cpart, cfNew, Rpart, Lpart, Rank, fail );
      else
        if cfEqu(4) > Dm
          CutIF.Lock = 'L+';
        elseif cfEqu(4) < Dm
          CutIF.Lock = 'R-';
        end
        Rpart = Rtemp;
        Lpart = Ltemp;
      end
    end
    
  end %IFBR 18
%% Simplify
  function [MST, IFT] = ReCutS( MST, IFT, Idx )
    MSL = MST.get(Idx);
    MSR = MST.get(Idx+1);
    Cobj = MST.get(MST.getparent(Idx));
    LIF = IFT.get(Idx-1);
    RIF = IFT.get(Idx+1);
    IF = IFT.get(Idx);
    if ~isempty(LIF)
      [Cobj, ~, ~] = Cobj.Cut(LIF.Equ);
    end
    if ~isempty(RIF)
      [~, Cobj, ~] = Cobj.Cut(RIF.Equ);
    end
    [Robj,Lobj,CutIF] = Cobj.Cut(IF.Equ,MSR,MSL,IF);
    MST = MST.set( Idx, Lobj );
    MST = MST.set( Idx+1, Robj );
    IFT = IFT.set( Idx, CutIF );
    MST = CutR(MST, IFT, Idx);
    MST = CutR(MST, IFT, Idx+1);
    
    function MST = CutR(MST, IFT, Idx)
      [MST, ~] = parallelCut(MST, IFT, Idx);
      cIdx = IFT.getchildren(Idx);
      n = numel(cIdx);
      for i = 1:n
        MST = CutR(MST, IFT, cIdx(i));
      end
    end
  end

  function [RVarVol, LVarVol] = moveVF(MST,IF,d)
    RVarVol = MST.VarVol; LVarVol = MST.VarVol;
    if IF.Rank
      IF = IF.subIF;
    end
    for j = 1:numel(IF)
      IFnow = IF(j);
      RVarVol = moveVarVol( RVarVol, IFnow, d);
      LVarVol = moveVarVol( LVarVol, IFnow, -d);
    end
  end
  
  function VarVol = moveVarVol( VarVol, IF, d )
    dV = IF.Area*d;
    VarVol(IF.RVar) = VarVol(IF.RVar) - dV;
    VarVol(IF.LVar) = VarVol(IF.LVar) + dV;
  end

  function [RwallEN, LwallEN] = movewall(IFT, IFs, bIFs, d)
    wallEN = IFT.EN;
    RwallEN = wallEN;
    LwallEN = wallEN;
    for i = 1:numel(bIFs)
      bIF = bIFs(i);
      mIFs = IFs.subnear(bIF);
      for j = 1:numel(mIFs)
        mIF = mIFs(j);
        sL = mIF.FxF(bIF);
        if ~isa(sL,'SegL')
          continue
        end
        rad = mIF.FvF(bIF);
        mIFdA = -d/tan(rad)*sL.long;
        bIFdA = -d/sin(rad)*sL.long;
        mIFdEN = mIFdA / mIF.Area * mIF.EN;
        bIFdEN = bIFdA / bIF.Area * bIF.EN;
        RwallEN = RwallEN + mIFdEN + bIFdEN;
        LwallEN = LwallEN - mIFdEN - bIFdEN;
      end
    end
  end
%% plot
  function data = replot(data, MST, IFT, EF, stress)
    if ~MSWeight
      delete(data.SC1);
      data.SC1 = MST.draw('Parent', ax.SC1);
      delete(data.IF1);
      data.IF1 = IFT.draw(IFtype, 'Parent', ax.IF1);
    else
      cidx = MST.getchildren(1);
      for i = 1:CN
        SCfield = ['SC',num2str(i)];
        IFfield = ['IF',num2str(i)];
        delete(data.(SCfield));
        delete(data.(IFfield));
        data.(SCfield) = MST.subtree(cidx(i)).draw('Parent',ax.(SCfield));
        axis(ax.(SCfield),axLim)
        data.(IFfield) = IFT.subtree(cidx(i)).draw(IFtype,'Parent',ax.(IFfield));
        axis(ax.(IFfield),axLim)
      end
    end
    
    E = EF(DirecE) / LoadDeg;
    P = MST.EffProp.PR(DirecE);
    eps = MST.EffProp.epsR(DirecE);
    VF = MST.VarVol / MST.Vol;
    G = getG(MST.VarVol,EF, stress);
    IFEN = IFT.EN;
    set(data.P, 'xdata', [data.P.XData, E], 'ydata', [data.P.YData, P]);
    set(data.eps, 'xdata', [data.eps.XData, E], 'ydata', [data.eps.YData, eps]);
    for i = 1:Const.vn
      set(data.V(i), 'xdata', [data.V(i).XData, t], 'ydata', [data.V(i).YData, VF(i)]);
    end
    set(data.G, 'xdata', [data.G.XData, t], 'ydata', [data.G.YData, G]);
    set(data.IFEN, 'xdata', [data.IFEN.XData, t], 'ydata', [data.IFEN.YData, IFEN]);
    set(data.sigma, 'string', ['sigma',dir,' = ',num2str(stress(DirecE)'/LoadDeg),' MPa']);
    set(data.EF,'string',['E',dir,' = ',num2str(EF(DirecE)/LoadDeg),' MV']);
    %% timer
    T = T + toc;
    hh = floor(rem(T,86400)/3600);
    mm = floor(rem(T,3600)/60);
    ss = floor(rem(T,60)/1);
    RT = toc*(tf-t)/dt;
    Rhh = floor(rem(RT,86400)/3600);
    Rmm = floor(rem(RT,3600)/60);
    Rss = floor(rem(RT,60)/1);
    set(hT, 'string', ['Time Passed: ',num2str(hh,'%02.0f'),':',...
      num2str(mm,'%02.0f'),':',num2str(ss,'%02.0f')] );
    set(hRT, 'string', ['Time Remaining: ',num2str(Rhh,'%02.0f'),':',...
      num2str(Rmm,'%02.0f'),':',num2str(Rss,'%02.0f')] );
  end

  function data = setfig(MST, IFT, EF, stress)
    n = ceil(sqrt(CN));
    SSize = 10;
    set(groot,'defaultAxesColorOrder',Const.varcolor);
    %% fcharts
    fcharts = figure('name','charts','numbertitle','off');
    ax.EN = subplot(2,2,1, 'nextplot', 'add', 'parent', fcharts);
    ax.P = subplot(2,2,3, 'nextplot', 'add', 'parent', fcharts);
    ax.eps = subplot(2,2,4, 'nextplot', 'add', 'parent', fcharts);
    ax.V = subplot(2,2,2, 'nextplot', 'add', 'parent', fcharts);
    FSize = 10;
    title(ax.EN, '$Energy-time$','interpreter','latex','fontsize',FSize);
    xlabel(ax.EN,'$time$','interpreter','latex','fontsize',FSize);
    ylabel(ax.EN,'$E(Jm^{-3})$','interpreter','latex','fontsize',FSize);
    title(ax.P, '$P-E$','interpreter','latex','fontsize',FSize);
    xlabel(ax.P,['$E_{',dir,'}(MVm^{-1})$'],'interpreter','latex','fontsize',FSize);
    ylabel(ax.P,['$D_{',dir,'}(Cm^{-2})$'],'interpreter','latex','fontsize',FSize);
    title(ax.eps, '$eps-E$','interpreter','latex','fontsize',FSize);
    xlabel(ax.eps,['$E_{',dir,'}(MVm^{-1})$'],'interpreter','latex','fontsize',FSize);
    ylabel(ax.eps,['$\varepsilon_{',dir,'}$(\%)'],'interpreter','latex','fontsize',FSize);
    title(ax.V,'$VF-time$','interpreter','latex','fontsize',FSize);
    xlabel(ax.V,'$time$','interpreter','latex','fontsize',FSize);
    ylabel(ax.V,'$Volume Fraction$','interpreter','latex','fontsize',FSize);
    %% SC,IF
    if ~MSWeight
      fSC = figure('Name','SCplot','NumberTitle','off');
      ax.SC1 = axes('parent', fSC, 'units','normalized',...
        'XLim',axLim(1:2),'YLim',axLim(3:4),'ZLim',axLim(5:6),...
        'xtick',[],'ytick',[],'ztick',[]);
      data.SC1 = MST.draw('Parent', ax.SC1);
      
      fIF = figure('Name', 'IFplot', 'numbertitle', 'off');
      ax.IF1 = axes('parent', fIF, 'units','normalized',...
        'XLim',axLim(1:2),'YLim',axLim(3:4),'ZLim',axLim(5:6),...
        'visible','off');
      data.IF1 = IFT.draw(IFtype, 'Parent', ax.IF1);
    else
      cidx = MST.getchildren(1);
      fSC = figure('Name','SCplot','NumberTitle','off');
      for i = 1:CN
        SCfield = ['SC',num2str(i)];
        ax.(SCfield) = subplot(n,n,i,'parent',fSC,'units','normalized');
        data.(SCfield) = MST.subtree(cidx(i)).draw('Parent',ax.(SCfield));
        axis(ax.(SCfield),axLim)
      end
      
      cidx = IFT.getchildren(1);
      fIF = figure('Name', 'IFplot', 'numbertitle', 'off');
      for i = 1:CN
        IFfield = ['IF',num2str(i)];
        ax.(IFfield) = subplot(n,n,i,'parent',fIF,'units','normalized');
        data.(IFfield) = IFT.subtree(cidx(i)).draw(IFtype,'Parent', ax.(IFfield));
        axis(ax.(IFfield),axLim)
      end
    end
    %% timers
    hT = uicontrol(fSC,'style','text','String','Time Passed: 00:00:00',...
      'FontWeight','bold','FontSize',10,'units','normalized','Position',[0.72 0.05 0.3 0.05]);
    hRT = uicontrol(fSC,'style','text','String','Time Remaining: 00:00:00',...
      'FontWeight','bold','FontSize',10,'units','normalized','Position',[0.7 0 0.3 0.05]);
    %% draw
    E = EF(DirecE) / LoadDeg;
    P = MST.EffProp.PR(DirecE);
    eps = MST.EffProp.epsR(DirecE);
    VF = MST.VarVol / MST.Vol;
    G = getG(MST.VarVol,EF, stress);
    IFEN = IFT.EN;
    data.P = plot(ax.P, E, P, 'b', 'marker','.','markersize',SSize);
    data.eps = plot(ax.eps, E, eps, 'b', 'marker','.','markersize',SSize);
    data.V = plot(ax.V, t, VF,'marker','.','markersize',SSize);
    data.G = plot(ax.EN, t, G, 'b', 'marker','.','markersize',SSize);
    data.IFEN = plot(ax.EN, t, IFEN, 'r', 'marker','.','markersize',SSize);
    data.sigma = uicontrol(fSC,'style','text','String',['sigma',dir,' = 0 MPa'],...
      'fontsize',10,'fontweight','bold','units','normalized','Position',[0 0 0.3 0.05]);
    data.EF = uicontrol(fSC,'style','text','String',['E',dir,' = 0 MV'],...
      'fontsize',10,'fontweight','bold','units','normalized','Position',[0 0.05 0.3 0.05]);
  end
end

% MST+IFT
% isinErr when the segL is very small, Tol?

% moveVF polyhedron Vol!!!
% moveIFwall,moveVF with bulk BIF!!!