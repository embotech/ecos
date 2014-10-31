% Driver for quadratic minimization with bounds and sparsity
function qpdriver(P, q, r, lb, ub, lambda, beta, alpha, print)

n = size(q, 1);
     
lambdaL1 = beta*lambda;
lambdaL2 = (1 - beta)*lambda;

rho = computeRho(P);

%L2 regularization

if lambdaL2 > 0
  P = P + lambdaL2*eye(n);
end

fprintf('variables %g lambdaL1 %g lambdaL2 %g beta %g alpha %g rho %g\n', n, lambdaL1, lambdaL2, beta, alpha, rho);

if lambdaL1 > 0
  Pt = [P zeros(n,n); zeros(n, 2*n)];
  qt = [q;lambdaL1*ones(n,1)];
  I  = eye(n);
  Aeq = [I -I;-I -I];
  beq = zeros(2*n, 1);
end

%Add MOSEK path
addpath /home/debasish/mosek/7/toolbox/r2013a/
tic;
if lambdaL1 > 0
   mosekx = quadprog(Pt, qt, Aeq, beq, [], [], [], []);
else
   mosekx = quadprog(P, q, [], [], [], [], lb, ub);
end
mosekTime = toc;
fprintf('mosek done\n');
mosekx = mosekx(1:n);

%Add ECOS path
addpath /home/debasish/ecos/matlab

if lambdaL1 > 0
  [ecosx, ~, ~, ~, ~, ecosTime] = ecosqp(Pt, qt, Aeq, beq, [], [], [], [], ecosoptimset('verbose', 0, 'feastol', 1e-8));
else 
  [ecosx, ~, ~, ~, ~, ecosTime] = ecosqp(P,q,[],[],[],[],lb,ub, ecosoptimset('verbose', 0, 'feastol', 1e-8));
end
fprintf('ecos done\n');
ecosx = ecosx(1:n);

tic;
[x history admmIters] = qpproximal(P, q, r, [], [], lb, ub, rho, alpha, lambdaL1);
admmTime = toc;

tic;
[xfast history accIters] = qpaccelerated(P, q, r, [], [], lb, ub, rho, lambdaL1);
admmFast = toc;

if print
fprintf('Mosek solution\n');
mosekx
fprintf('Proximal solution\n');
x
end

fprintf('mosek-ecos norm %g\n', norm(mosekx - ecosx, Inf));
fprintf('mosek-admm norm %g\n', norm(mosekx - x, Inf));
fprintf('mosek-admmfast norm %g\n', norm(mosekx - xfast, Inf));
fprintf('mosek %g ecos %g admm %g iters %g admmfast %g accIters %g\n', mosekTime, ecosTime, admmTime, admmIters, admmFast, accIters);

end
