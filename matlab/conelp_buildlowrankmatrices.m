function [S,V] = conelp_buildlowrankmatrices(scalings,n,p,dims)
% Returns the low rank part coming from second-order cone scalings.


k = length(dims.q);
mtilde = dims.l + sum(dims.q) + k;

V = zeros(n+p+mtilde,k);
s = NaN(k,1);

for j = 1:k
    V(n+p+dims.l+sum(dims.q(1:j-1))+j,j) = 1;
    s(j) = scalings.q(j).beta;
end

S = diag(s);