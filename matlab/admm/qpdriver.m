% Standard form QP example
function qpdriver(n, m)
fprintf('variables %g equality %g\n', n, m)

randn('state', 0);
rand('state', 0);

% generate a well-conditioned positive definite matrix
% (for faster convergence)
P = rand(n, n);
P = P'*P;
P = (P + P')/2;
%P = P + P';
%[V D] = eig(P);
%P = V*diag(1+rand(n,1))*V';

q = randn(n,1);
r = randn(1);

l = randn(n,1);
u = randn(n,1);
lb = min(l,u);
ub = max(l,u);

%Add MOSEK path
addpath /home/debasish/mosek/7/toolbox/r2013a/
tic;
mosekx = quadprog(P, q, [], [], [], [], lb, ub);
mosekTime = toc;
fprintf('mosek done\n');

%Add ECOS path
addpath /home/debasish/ecos/matlab
[ecosx, ~, ~, ~, ~, ecosTime] = ecosqp(P,q,[],[],[],[],lb,ub, ecosoptimset('verbose', 0, 'feastol', 1e-8));

tic;
[x history] = qpproximal(P, q, r, lb, ub, 1.0, 1.0);
admmTime = toc;

tic;
[xfast history] = qpaccelerated(P,q,r,lb,ub,1.0,1.0);
admmFast = toc;

fprintf('mosek-ecos norm %g\n', norm(mosekx - ecosx));
fprintf('mosek-admm norm %g\n', norm(mosekx - x));
fprintf('mosek-admmfast norm %g\n', norm(mosekx - xfast));
fprintf('mosek %g ecos %g admm %g admmfast %g\n', mosekTime, ecosTime, admmTime, admmFast);
