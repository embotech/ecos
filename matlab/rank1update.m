function [Lbar,Dbar] = rank1update(L,D,alpha,z)
% Updates Cholesky factors of A + alpha*z*z' given A = LDL'.
% Implements C1 method from the paper "Method for modifying matrix
% factorizations" by Gill, Golub, Murray and Saunders, p. 516.

n = size(L,1);

w = z;
d = diag(D);
dbar = d;
Lbar = L;
for j = 1:n
    if( w(j) ~= 0 )
        p = w(j);
        dbar(j) = d(j) + alpha*p^2;
        beta = p*alpha / dbar(j);
        alpha = d(j)*alpha / dbar(j);
        w(j+1:end,1) = w(j+1:end,1) - p*L(j+1:end,j);
        Lbar(j+1:end,j) = L(j+1:end,j) + beta*w(j+1:end,1);
    end
end
Dbar = diag(dbar);
