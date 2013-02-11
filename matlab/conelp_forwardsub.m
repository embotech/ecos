function x = conelp_forwardsub(L,b)
% Solves Lx = b by forward substitution

n = size(L,1);
assert(size(L,2) == n,'L must be square');
assert(size(b,1) == n,'dimension mismatch');
assert(size(b,2) == 1,'b must be a column vector');

x = NaN(n,1);
for i = 1:n
    a = 0;
    for j = 1:i-1
        a = a + L(i,j)*x(j);
    end
    x(i) = (b(i) - a) / L(i,i);
end