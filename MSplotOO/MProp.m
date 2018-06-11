classdef MProp
    %% MPROP  Material properties generator for each variant
    
    properties (Hidden, Constant)
        Ps = 0.27;
        epss = 0.8/100;
        refP = getrefP('t');
    end
    
    properties (SetAccess = protected)
        PR = zeros(1,3);
        epsR = zeros(1,6);
        k = zeros(3,3);
        s = zeros(6,6);
        d = zeros(3,6);
        c = zeros(6,6);
        e = zeros(3,6);
    end
    
    methods
        function obj = MProp(var)
            if nargin
                assert( isnumeric(var) && all(~mod(var,1)), 'Invalid input for Mprop.');
                n = length(var);
                obj(n) = MProp;
                for i = 1:n
                    obj(i) = prepare(obj(i),var(i));
                end
            end            
        end
        
        function objA = plus( objA, objB )
            %% 
            
            if isa( objA, 'MProp' ) && isa( objB, 'MProp' )
                objA.PR = objA.PR + objB.PR;
                objA.epsR = objA.epsR + objB.epsR;
                objA.k = objA.k + objB.k;
                objA.s = objA.s + objB.s;
                objA.d = objA.d + objB.d;
                objA.c = objA.c + objB.c;
                objA.e = objA.e + objB.e;
            end    
            
        end
        
        function obj = mtimes( obj, num )
            %%
            
            if isa( obj, 'MProp' ) && isscalar( num )
                obj.PR = obj.PR*num;
                obj.epsR = obj.epsR*num;
                obj.k = obj.k*num;
                obj.s = obj.s*num;
                obj.d = obj.d*num;
                obj.c = obj.c*num;
                obj.e = obj.e*num;
            end
            
        end
    end
    
end

function refP = getrefP( type )
%make refP=the array of variant property
refP = [];
switch type
    case 'r' % Rhombohedral polar Direction.
        Pt=1;   vn=8;% Number of variants
        P=[ 1,   1,  1;
            -1,   1,  1;
            -1,  -1,  1;
            1,  -1,  1;  ]*Pt;% Treat '+' and '-' direction as a pair.
        n=0;    d=1;    et=1;
        e=[ n,   d,  d,  n,  d,  n;
            n,  -d, -d,  n,  d,  n;
            n,   d, -d,  n, -d,  n;
            n,  -d,  d,  n, -d,  n;]*et;
    case 't' % Tetragonal polar direction.
        Pt=1;   vn=6;% Number of variants
        P=[ 1,  0,  0;
            0,  1,  0;
            0,  0,  1]*Pt;
        alpha=-0.5; beta=1; et=1;
        e=[ beta,    0,  0,  alpha,  0,  alpha;
            alpha,   0,  0,  beta,   0,  alpha;
            alpha,   0,  0,  alpha,  0,  beta;  ]*et;
end
for id = 1: vn/2 % Put P and e together.
    refP = cat(1,refP, [P(id,:),e(id,:)]);
    refP = cat(1,refP, [-P(id,:),e(id,:)]);
end
end

function obj = prepare(obj,V)
% Set properties of BaTiO3

%% Material constants
c11=275e9;   c12=179e9;   c13=152e9; %original
c33=165e9;   c44=54.3e9;  c66=113e9;
e15 =21.3;   e31=-2.69;   e33=3.65;
kapa0=8.854e-12;
kapa11=1970*kapa0;  kapa33=109*kapa0;
beta=-0.2222;
%% Properties of each crystal vsriant
switch V
    case { 1, 2 }
        c=[c33 c13 c13 0   0   0;
            c13 c11 c12 0   0   0;
            c13 c12 c11 0   0   0;
            0   0   0   c66 0   0;
            0   0   0   0   c44 0;
            0   0   0   0   0   c44];
        k=[kapa33 0 0; 0 kapa11 0; 0 0 kapa11];
        e=[e33 e31 e31 0   0   0;
            0   0   0   0   0   e15;
            0   0   0   0   e15 0 ]*(-1)^mod(V-1,2);
        PR = [1;0;0]*(-1)^mod(V-1,2)*obj.Ps;
        epsR = [1;beta;beta;0;0;0]*obj.epss;
    case { 3, 4 }
        c=[c11 c13 c12 0   0   0;
            c13 c33 c13 0   0   0;
            c12 c13 c11 0   0   0;
            0   0   0   c44 0   0;
            0   0   0   0   c66 0;
            0   0   0   0   0   c44];
        k=[kapa11 0 0; 0 kapa33 0; 0 0 kapa11];
        e=[0   0   0   0   0   e15;
            e31 e33 e31 0   0   0  ;
            0   0   0   e15 0   0 ]*(-1)^mod(V-1,2);
        PR = [0;1;0]*(-1)^mod(V-1,2)*obj.Ps;
        epsR = [beta;1;beta;0;0;0]*obj.epss;
    case { 5, 6 }
        c=[c11 c12 c13 0   0   0;
            c12 c11 c13 0   0   0;
            c13 c13 c33 0   0   0;
            0   0   0   c44 0   0;
            0   0   0   0   c44 0;
            0   0   0   0   0   c66];
        k=[kapa11 0 0; 0 kapa11 0; 0 0 kapa33];
        e=[0   0   0   0   e15 0;
            0   0   0   e15 0   0;
            e31 e31 e33 0   0   0]*(-1)^mod(V-1,2);
        PR = [0;0;1]*(-1)^mod(V-1,2)*obj.Ps;
        epsR = [beta;beta;1;0;0;0]*obj.epss;
    otherwise
        error( 'Invalid "V".' );
end
s = inv(c);
d = e/c;
%% Material properties transformation
refP = obj.refP;
refP(:,4:9) = [ refP(:,4), refP(:,7), refP(:,9), refP(:,8), refP(:,6), refP(:,5) ]; %temporary modify
refP = abs(refP');
obj.PR = ( PR.*refP(1:3,V) )';
obj.epsR = ( epsR.*refP(4:9,V) )';
obj.k = k;
obj.s = s;
obj.d = d;
obj.c = c;
obj.e = e;

end