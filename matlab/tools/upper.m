function U = u(A)
% returns upper triangular part of symmetric matrix A

n = size(A,1);
assert(size(A,2)==n,'A must be square');

U = zeros(n);

for i = 1:n
    U(1:i,i) = A(1:i,i);
end