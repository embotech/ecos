function x = conelp_backwardsub(U,b)
% Solves Ux = b by backward substitution

n = size(U,1);
assert(size(U,2) == n,'L must be square');
assert(size(b,1) == n,'dimension mismatch');
assert(size(b,2) == 1,'b must be a column vector');

x = NaN(n,1);
for i = n:-1:1
    a = 0;
    for j = n:-1:i+1
        a = a + U(i,j)*x(j);
    end
    x(i) = (b(i) - a) / U(i,i);
end