function x = myldl(A,alpha,z,b)
% LDL solver using rank-one modifications.
%
% Solves the system (A + alpha*zz')x = b by LDL factorization
% of A and low-rank forward/backward substitutions.

n = size(A,2);
assert(length(z) == n,'Dimension mismatch v and A');
assert(size(b,1) == size(A,1), 'Dimension mismatch A and b');

% get an LDL of A and de-sparsify it
[L,D] = ldlsparse(sparse(A),1:n);
L = full(L) + eye(n);
D = full(D);

% rank-1-update Cholesky factors
[Lbar,Dbar] = rank1update(L,D,alpha,z);
%[Ltrue,Dtrue] = ldlsparse(sparse(A+alpha*v*v'),1:n);
%Ltrue = Ltrue + eye(n);
%norm(Ltrue-Lbar)
%norm(Dbar-Dtrue)

% solve
%y = conelp_forwardsub_bydiag_r1(L,D,alpha,z,b);
%y = conelp_byDiag(Dbar,z);
%x = conelp_backwardsub(Lbar',y);

[x,y,z] = conelp_forwardsub_bydiag_backwardsub_r1(L,D,alpha,z,b);
%y = conelp_byDiag(Dbar,z);
x = conelp_backwardsub(Lbar',y);