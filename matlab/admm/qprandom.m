%Generate random QPs with bound and sparsity
function qprandom(n, lambda, beta, alpha)

randn('state', 10);
rand('state', 10);

P = rand(n, n);
P = P'*P;

q = randn(n,1);
r = randn(1);

l = randn(n,1);
u = randn(n,1);

lb = min(l,u);
ub = max(l,u);

qpdriver(P, q, r, lb, ub, lambda, beta, alpha, false);
end
