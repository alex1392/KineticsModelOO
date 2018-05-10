function [STree,Badness,numAnswer] =BuildTree(obj)
Var = obj.Var;
DOF = obj.DOF;
Tol = obj.Tol;
if isempty(DOF)
  DOF = 1;
end
FavoredP=[];
for idx = 0 : length(Var)-1 %make the bottom node
  fracIdx=de2x(idx,2,length(DOF));
  fraction=1;
  for jdx = 1 : length(fracIdx)
    fraction=fraction*(DOF(jdx)*(-2*fracIdx(jdx)+1)+fracIdx(jdx)); %calculate volume fraction
  end
  MP = MProp;
  FavoredP = cat(2 , FavoredP , Node(fraction, MP.refP(Var(idx+1),:)));
end
[STree,Badness,numAnswer] = Construct(FavoredP);  %Construct the upper node of the tree
if iscell( STree )% Generalize the Tree
  STree = STree{randi(length(STree))};
  %STree = STree{3}; %for [1 2 3 4]
end
  
if isnan( STree(1).IfNormal(1) )
  STree = InverseArray( STree );
end

%%
  function [Answer,Badness,numAnswer] = Construct(inputFP)
    
    FP=inputFP; %input Nodes
    [Answer,NeedModifP,DoubleStrain,Success] = Laminate(FP); % Laminate the bottom layer into a tree
    Badness=1e6;
    numAnswer=1;
    % NeedModifP = 0 indicates there is no 180 degree domain wall.
    % DoubleStrain = 0 indicates there is no domain wall has two possible
    % interface normals selected by strain method.
    if Success
      Answer = InverseArray(Answer); % Inverse array for better data structure.
      if DoubleStrain==0 && NeedModifP==0
        Badness = BadnessSolver(Answer);
      else
        [Answer,Badness,numAnswer] = OptIfnormal(Answer,DoubleStrain,NeedModifP);
      end
      if numAnswer>1
        for Idx1=1:numAnswer
          Answer{Idx1} = InverseArray(Answer{Idx1});
        end
      else
        Answer = InverseArray(Answer); % Inverse array for better data structure.
      end
    end
    
    % Construct Inner Function
    function [Answer,Badness,numAnswer] = OptIfnormal(Answer,DoubleStrain,NeedModifP)
      numAnswer=0;
      ABest = Answer;
      Badness=1e6;
      % Badness =  BadnessSolver(Answer);
      if DoubleStrain>0
        DoubleStrain = 0;
        index=zeros(1,size(Answer,2));
        for i = 1:size(Answer,2)
          if strcmp(Answer(i).Method,'Stran2')
            DoubleStrain = DoubleStrain+1;
            index(DoubleStrain)=i;
          end
        end
        Astandard=Answer;
        Answer=[];
        for i =1 : 2^DoubleStrain
          Atemp =Astandard;
          pattern = zeros(1,DoubleStrain);
          n = de2x(2^DoubleStrain-i,2,DoubleStrain);
          pattern(1,1:size(n,2))=n;
          for j = 1:DoubleStrain
            Atemp(index(j)).IfNormal(:,1)=Atemp(index(j)).IfNormal(:,pattern(j)+1);
          end
          if NeedModifP>0
            [Atemp,numAnswer] = ModifyPP(Atemp,numAnswer);
          end
          if numAnswer==1
            Btemp =  BadnessSolver(Atemp);
            if Btemp<Badness || (abs(Btemp-Badness)<Tol)
              ABest = Atemp;
              Badness =  Btemp;
            end
            if abs(Badness)<Tol
              if isempty(Answer)
                Answer= ABest;
              else
                Answer=cat(2,Answer,ABest);
              end
            end
          elseif numAnswer>1
            if isempty(Answer)
              Answer=Atemp;
            else
              Answer=cat(2,Answer,Atemp);
            end
            Badness=0;
          end
        end
        if numAnswer==1 && isempty(Answer)
          Answer = ABest;
        end
      else
        [Atemp,numAnswer] = ModifyPP(Answer,numAnswer);
        if numAnswer ==1
          Btemp =  BadnessSolver(Atemp);
          if Btemp<Badness
            Answer = Atemp;
            Badness =  Btemp;
            return
          end
        else
          Answer = Atemp;
          Badness=0;
        end
      end
    end
    
    function  [Result,NeedModifP,DoubleStrain,Success]= Laminate(NodeSet)
      global FixedInf FixNodeNum
      NeedModifP = 0; % 1 when there is any 180 degree domain wall happened
      DoubleStrain = 0; % 1 when there is any laminate has parelle polarization but different strain state.
      Success=1;
      % Check if the bottom layer is valid.
      Size = size(NodeSet,2);
      Valid = 0;
      for n = 0 : 7
        if Size == 2^n
          Valid = 1;
          break;
        end
      end
      for i=1:size(NodeSet,2)
        if NodeSet(i).Lambda<=0
          Valid = false;
        end
      end
      
      if ~Valid % if not valid, return
        Result(1).IfNormal = [];
        Success=0;
        return;
      else % if it is valid, laminate them.
        start=1;
        for next = log2(Size)-1: -1: 0
          LowLayer = NodeSet( start: size(NodeSet,2));
          for n = 1 : 2^next
            % Calculate all the propertise in a volume average sense.
            upLayer(n).Lambda = LowLayer(2*n-1).Lambda+ LowLayer(2*n).Lambda; % Volume fraction
            upLayer(n).P      = (LowLayer(2*n-1).Lambda*LowLayer(2*n-1).P+ LowLayer(2*n).Lambda*LowLayer(2*n).P) / upLayer(n).Lambda; % Polarization
            upLayer(n).e      = (LowLayer(2*n-1).Lambda*LowLayer(2*n-1).e+ LowLayer(2*n).Lambda*LowLayer(2*n).e) / upLayer(n).Lambda; % Strain
            % Calculate interface normal vector.
            [normal ,text,NeedModifP,DoubleStrain] = InterfaceNormal(LowLayer(2*n-1:2*n),NeedModifP,DoubleStrain);
            upLayer(n).IfNormal = normal; % Interface normal
            upLayer(n).Method = text; % Method of getting interface normal
            if isequal(upLayer(n).IfNormal(:,1),[Inf;Inf;Inf])&& ...
                ~strcmp(upLayer(n).Method,'P==P  ') && ~strcmp(upLayer(n).Method,'P//P  ') &&...
                ~strcmp(upLayer(n).Method,'Polar ') % If any interface normal does not exist, return.
              NodeSet(1).Lambda = -99;
              Result = NodeSet;
              Success = 0;
              return
            end
          end
          NodeSet = cat(2, NodeSet , upLayer(1:2^next));
          start = start + 2^(next+1);
        end
        if ~isempty(FixedInf)
          NodeSet(FixNodeNum).IfNormal = FixedInf; % Interface normal
          NodeSet(FixNodeNum).Method = 'Fixed '; % Method of getting interface normal
        end
        Result = NodeSet;
      end
    end
    
    function [IfNormal,text,NeedModifP,DoubleStrain] = InterfaceNormal(NodeSet,NeedModifP,DoubleStrain)
      PDifference=NodeSet(1).P-NodeSet(2).P;
      eDifference=NodeSet(1).e-NodeSet(2).e;
      IfNormal(:,:)=[[Inf;Inf;Inf],[Inf;Inf;Inf]];
      isZeroM = 0;
      text = 'Fail ';
      
      for i=1:3
        for j = 1:3
          isZeroM = isZeroM+ abs(eDifference(i,j));
        end
      end
      
      if isZeroM< Tol
        % 1. When P = P
        if abs(dot(NodeSet(1).P/norm(NodeSet(1).P),NodeSet(2).P/norm(NodeSet(2).P)) -1) < Tol...
            || abs(norm(NodeSet(1).P)+norm(NodeSet(2).P))<Tol
          if all(isnan(NodeSet(1).IfNormal(:,1))==1) && all(isnan(NodeSet(2).IfNormal(:,1))==1)
            text = 'P==P  ';
            IfNormal(:,1) = [0.0;0.0;0.0];
          elseif norm(NodeSet(1).P-NodeSet(2).P)<Tol
            text = 'P==P  ';
            IfNormal(:,1) = [0.0;0.0;0.0];
          else
            text = 'P//P  ';
            IfNormal(:,1) = NodeSet(1).P/norm(NodeSet(1).P);
          end
          NeedModifP = NeedModifP+1;
          % 2. When P = -P
        elseif abs(dot(NodeSet(1).P/norm(NodeSet(1).P),NodeSet(2).P/norm(NodeSet(2).P)) +1) < Tol
          if norm(cross(NodeSet(1).P, NodeSet(2).P ))<Tol
            text = 'P//P  ';
            IfNormal(:,1) = NodeSet(1).P/norm(NodeSet(1).P);
            NeedModifP = NeedModifP+1;
          end
          % 3. Make (P1-P2). n == 0 .
        else
          PlaneNormal = NodeSet(1).P-NodeSet(2).P;
          text ='Polar ';
          NeedModifP = NeedModifP+1;
          IfNormal(:,1)=PlaneNormal/norm(PlaneNormal);
        end
      else
        % 4. n = +- sqrt(-eValue1 / eValue2)* e1 + e2
        if abs(det(eDifference)) < Tol && (trace(eDifference)^2-trace(eDifference^2))< -Tol
          D=eig(eDifference);
          [V,~]=eig(eDifference);
          % Assign the negative eigenvalue to eValue1, and the positive one is eValue2.
          for n=1:3
            if D(n,1) < -Tol
              eValue1=D(n,1);
              eVector1=V(:,n);
            elseif D(n,1) > Tol
              eValue2=D(n,1);
              eVector2=V(:,n);
            end
          end
          n1=(   sqrt(-1*eValue1/eValue2)* eVector1 + eVector2);
          n2=(-1*sqrt(-1*eValue1/eValue2)* eVector1 + eVector2);
          % Check if the normal is compatible for polarization.
          if abs(dot(PDifference,n1)) < Tol && abs(dot(PDifference,n2)) > Tol % n1 fits both criteria of compatibility
            IfNormal(:,1)=n1;
          elseif abs(dot(PDifference,n1)) > Tol && abs(dot(PDifference,n2)) < Tol % n2 fits both criteria of compatibility
            IfNormal(:,1)=n2;
          elseif abs(dot(PDifference,n1)) < Tol && abs(dot(PDifference,n2)) < Tol % n1 and n2 fit both criteria of compatibility
            IfNormal=[n1,n2];
          else % Nothing fits criteria of compatibility
            IfNormal(:,1)=[Inf;Inf;Inf];
          end
          
          % Check again by the function of CheckStrain()
          if ~isequal(IfNormal(:,1),[Inf;Inf;Inf]) && norm(IfNormal(:,1))>Tol && isequal(IfNormal(:,2),[Inf;Inf;Inf])  % If only one normal vector is chosen.
            Compatible=CheckStrain(IfNormal(:,1),eDifference);
            if ~Compatible
              IfNormal(:,1)=[Inf;Inf;Inf];
              % It sould not be happended....
            else
              text = 'Strain';
            end
          elseif ~isequal(IfNormal(:,1),[Inf;Inf;Inf]) && norm(IfNormal(:,1))>Tol && ...
              ~isequal(IfNormal(:,2),[Inf;Inf;Inf]) && norm(IfNormal(:,2))>Tol  % If a pair normal vectors are chosen.
            Compatible=CheckStrain(IfNormal(:,1),eDifference);
            if ~Compatible
              IfNormal(:,1) = IfNormal(:,2);
              IfNormal(:,2) = [Inf;Inf;Inf];
              Compatible=CheckStrain(IfNormal(:,1),eDifference);
              if ~Compatible % It sould not be happended....
                IfNormal(:,1) = [Inf;Inf;Inf];
              else
                text = 'Strain';
              end
            else
              Compatible=CheckStrain(IfNormal(:,2),eDifference);
              if ~Compatible
                IfNormal(:,2) = [Inf;Inf;Inf];
                text = 'Strain';
              else
                text = 'Stran2';
                DoubleStrain = DoubleStrain+1;
              end
            end
          end
        else
          IfNormal(:,1)=[Inf;Inf;Inf];
        end
        i = 0;
        for n = 1:2
          if ~isequal(IfNormal(:,n),[Inf;Inf;Inf])&& norm(IfNormal(:,n))>Tol
            IfNormal(:,n)=IfNormal(:,n)/norm(IfNormal(:,n));
            i = i+1;
            if acosd(dot(IfNormal(:,n),[0;0;1]))>acosd(dot(-IfNormal(:,n),[0;0;1]))
              IfNormal(:,n) = IfNormal(:,n)*(-1);
            end
          end
        end
      end
      
      % If the interface normal still doesn't exist, construction is fail.
      if strcmp(text,'Fail  ')
        IfNormal(1,:)=[Inf;Inf;Inf];
        IfNormal(2,:)=[Inf;Inf;Inf];
      end
    end
    
    function [ABest,numAnswer] = ModifyPP(Answer,numAnswer)
      NeedModifP = 0;
      index = zeros(1,size(Answer,2));
      for i = 1: size(Answer,2)
        if strcmp(Answer(i).Method,'P==P  ') || strcmp(Answer(i).Method,'P//P  ') || strcmp(Answer(i).Method,'Polar ')
          NeedModifP = NeedModifP+1;
          index(NeedModifP) = i;
        end
      end
      
      rP(:,1)=[-1.0/sqrt(2); 1.0/sqrt(2); 0.0];
      rP(:,2)=[ 1.0/sqrt(2); 1.0/sqrt(2); 0.0];
      rP(:,3)=[ 0.0;-1.0/sqrt(2); 1.0/sqrt(2)];
      rP(:,4)=[ 0.0; 1.0/sqrt(2); 1.0/sqrt(2)];
      rP(:,5)=[-1.0/sqrt(2); 0.0; 1.0/sqrt(2)];
      rP(:,6)=[ 1.0/sqrt(2); 0.0; 1.0/sqrt(2)];
      rP(:,7)=[ 1.0; 0.0; 0.0];
      rP(:,8)=[ 0.0; 1.0; 0.0];
      rP(:,9)=[ 0.0; 0.0; 1.0];
      
      ABest = cell(1, 9^NeedModifP);
      numAnswerNow=0;
      for i = 1 : 9^NeedModifP
        TryThis = 1;
        Atemp = Answer;
        pattern = de2x(9^NeedModifP-i , 9 , NeedModifP);
        for j = 1: NeedModifP
          if abs(dot(rP(:,pattern(j)+1),Answer(index(j)).IfNormal(:,1))) < Tol
            Atemp(index(j)).IfNormal(:,1) = rP(:,pattern(j)+1);
          else
            TryThis=0;
            break
          end
        end
        if TryThis
          Btemp =  BadnessSolver(Atemp);
          if abs(Btemp) < Tol
            numAnswerNow = numAnswerNow+1;
            ABest{numAnswerNow} = Atemp;
          end
        end
      end
      ABest(numAnswerNow+1 : end) = [];
      
      if numAnswerNow==1
        ABest=ABest{1};
      elseif isempty(ABest)
        ABest = Answer;
      end
      numAnswer=numAnswer+numAnswerNow;
    end
    
    function Compatible = CheckStrain(n,M)
      
      A=[ 2*n(1),   0,        0;
        n(2),     n(1),     0;
        n(3),     0,        n(1);
        0,        2*n(2),   0;
        0,        n(3),     n(2);
        0,        0,        2*n(3)]*0.5;
      
      b=[M(1,1);M(1,2);M(1,3);M(2,2);M(2,3);M(3,3)];
      %         a= SolveAxb(A,b);
      a= A\b;
      
      % if a exist, the normal is stress-compatible.
      if isempty(a)
        Compatible = 0;
      else
        Compatible = 1;
      end
      
    end %CheckStrain 14(35)
    
    function Badness =  BadnessSolver(ns)
      
      Badness =0;
      sumNode=0;
      for n = 0 : 7
        sumNode = sumNode+2^n;
        if size(ns,2) == sumNode
          break;
        end
      end
      HRank = n; % Overall rank obtained.
      
      if ns(length(ns)).Lambda == -99
        Badness = 10^6;
      else
        for i=1:2^(HRank-1)-1
          n3 = ns(i).IfNormal(:,1);
          Rank = floor(log(i*1.0)/log(2.0));
          for r = 1: HRank -Rank-1
            for j = 1 : 2^(r-1)
              % Ranking the orientation
              n1 = ns(i*2^r+j-1).IfNormal(:,1);
              n2 = ns(i*2^r+j-1+2^(r-1)).IfNormal(:,1);
              n1 = n1 - n3*dot(n1,n3);
              n2 = n2 - n3*dot(n2,n3);
              if norm(n1)>Tol && norm(n2)>Tol
                n1 = n1 / norm(n1);
                n2 = n2 / norm(n2);
                Theta1 = dot(n1,n2);
                if abs(abs(Theta1)-1)<Tol
                  Theta1 = Mysign(Theta1);
                end
                Theta1 = acos(Theta1);
                Theta2 = dot(-n1,n2);
                if abs(abs(Theta2)-1)<Tol
                  Theta2 = Mysign(Theta2);
                end
                Theta2 = acos(Theta2);
                if Theta1 > Theta2
                  Theta1 = Theta2;
                end
                MisAlighment=abs(dot(ns((i*2^r+j-1)*2).P-ns((i*2^r+j-1+2^(r-1))*2).P,n3));
                MisAlighment=MisAlighment+abs(dot( ns((i*2^r+j-1)*2+1).P-ns((i*2^r+j-1+2^(r-1))*2+1).P,n3));
                Compatible = CheckStrain(n3,ns((i*2^r+j-1)*2).e-ns((i*2^r+j-1+2^(r-1))*2).e);
                if Compatible % Penalize incompatible domain pairs
                  StrainC=0;
                else
                  StrainC=10;
                end
                Compatible = CheckStrain(n3,ns((i*2^r+j-1)*2+1).e-ns((i*2^r+j-1+2^(r-1))*2+1).e);
                if Compatible
                  StrainC=StrainC+0;
                else
                  StrainC=StrainC+10;
                end
                PartBadness = abs((Theta1/(pi/2)+...
                  abs(ns((i*2^r+j-1)*2).Lambda/ns((i*2^r+j-1)*2+1).Lambda-(ns((i*2^r+j-1+2^(r-1))*2).Lambda/ns((i*2^r+j-1+2^(r-1))*2+1).Lambda)) +...
                  abs(MisAlighment)+abs(StrainC))* 5^(floor(log(i*1.0)/log(2.0)))* ns(i).Lambda);
              else
                PartBadness = 0;
              end
              Badness = Badness + PartBadness;
            end
          end
        end
      end
      
      % BadnessSolver Inner Function
      function result=Mysign(x)
        
        if x<-Tol
          result = -1;
        elseif x>Tol
          result = 1;
        else
          result = 0;
        end
        
      end %Mysign 6
      
    end %BadnessSolver 63(6)
    
  end

  function N = Node(Lambda , Variant)
    %make tree nodes for a variant
    N.Lambda = Lambda;                              % Volume fraction
    N.P = (Variant(1,1:3))';                        % Polarization
    N.e = [Variant(1,4),Variant(1,5),Variant(1,6);  % Strain
      Variant(1,5),Variant(1,7),Variant(1,8);
      Variant(1,6),Variant(1,8),Variant(1,9)];
    N.IfNormal = nan(3,1);                     % Default interface normal
    N.Method = '---';                               % Default method of getting interface normal
  end

  function s=de2x(n,x,L)
    %x-nary commutation  n=number, x=x-nary,L=zero fill to L digit
    
    s=[];
    while n>1
      s=cat(2,mod(n,x),s);
      n=floor(n/x);
    end
    if n==1
      s=[n s];
    end
    if isempty(s)
      s=0;
    end
    s=InverseArray(s);
    s(length(s)+1:L)=0;
  end

  function y = InverseArray(x)
    
    n = length(x);
    y = x;
    for i = 1: n
      y(i) = x(n - i+1);
    end
    
  end

end

