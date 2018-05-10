function [A, B] = swap(A, B)
temp = A;
A = B;
B = temp;
end