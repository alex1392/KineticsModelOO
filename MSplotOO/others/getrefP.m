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