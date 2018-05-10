function dis = E_E(EA,EB)
assert(isparallel(EA(1:3),EB(1:3)),'The Planes should be parallel.');
n = size(EB,1);
dis = zeros(n,1);
for i = 1:n
    dis(i) = abs(EA(4) - EB(4));
end
end