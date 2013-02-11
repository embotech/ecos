function [x,y,z] = conelp_forwardsub_bydiag_backwardsub_r1(L,D,alpha,v,b)
% Solves DLx = b by forward substitution with rank 1 update.
% L is the LDL factor of A, and we want to solve A + alpha*v*v'.

n = size(L,1);
assert(size(L,2) == n,'L must be square');
assert(size(b,1) == n,'dimension mismatch');
assert(size(b,2) == 1,'b must be a column vector');

x = NaN(n,1);

w = v;
d = diag(D);
dbar = NaN(n,1);
beta = NaN(n,1);
alpha = [alpha; NaN(n-1,1)];
p = NaN(n,1);
y = NaN(n,1);
z = NaN(n,1);
x = NaN(n,1);
for j = 1:n   
    p(j) = w(j);
    dbar(j) = d(j) + alpha(j)*p(j)^2;
    beta(j) = p(j)*alpha(j) / dbar(j);
    if( j<n), alpha(j+1) = d(j)*alpha(j) / dbar(j); end
    w(j+1:end,1) = w(j+1:end,1) - p(j)*L(j+1:end,j);  
    
    z(j) = b(j);
    y(j) = z(j) / dbar(j);
    b(j+1:end) = b(j+1:end) - z(j).*(L(j+1:end,j) + beta(j)*w(j+1:end,1));
end

yy = y;
x=y;
% for j = n:-1:1
%     x(j) = yy(j);
%     w(n:-1:j+1) = w(n:-1:j+1) + p(j)*L(n:-1:j+1,j);
%     yy(j-1:n:-1:j+1) = yy(j-1:-1:1) - x(j).*(L(j-1:-1:1,j) + beta(j)*w(j-1:-1:1,1));
% end

