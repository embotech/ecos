function x = conelp_forwardsub_r1(L,D,alpha,z,b)
% Solves Lx = b by forward substitution with rank 1 update.
% L is the LDL factor of A, and we want to solve A + alpha*v*v'.

n = size(L,1);
assert(size(L,2) == n,'L must be square');
assert(size(b,1) == n,'dimension mismatch');
assert(size(b,2) == 1,'b must be a column vector');

x = NaN(n,1);

w = z;
d = diag(D);

for j = 1:n   
    
    p = w(j);
    dbar(j) = d(j) + alpha*p^2;
    beta = p*alpha / dbar(j);
    alpha = d(j)*alpha / dbar(j);
    w(j+1:end,1) = w(j+1:end,1) - p*L(j+1:end,j);  
    
    x(j) = b(j);
    b(j+1:end) = b(j+1:end) - x(j).*(L(j+1:end,j) + beta*w(j+1:end,1));

end
