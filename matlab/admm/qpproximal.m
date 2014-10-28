function [z, history, admmIters] = qpproximal(P, q, r, Aeq, beq, lb, ub, rho, alpha, lambdaSparse)

% quadprog  Solve standard form box-constrained QP via ADMM
%
% [x, history] = quadprog(P, q, r, lb, ub, rho, alpha)
%
% Solves the following problem via ADMM:
%
%   minimize     (1/2)*x'*P*x + q'*x + r
%   subject to  Aeqx <= beq
%               lb <= x <= ub
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
    
QUIET    = 1;

MAX_ITER = max(400, 20*size(P,1));
ABSTOL   = 1e-8;
RELTOL   = 1e-4;

n = size(P,1);
equalities = size(Aeq, 1);

x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

admmIters = 0;

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

if equalities > 0
   H = [P + rho*eye(n), Aeq'; Aeq, zeros(equalities)];
   [L, U] = lu(H);
else
   R = chol(P + rho*eye(n));
end

for k = 1:MAX_ITER
    if equalities > 0
        scale = [rho*(z - u) - q;beq];
        xlambda = U \ (L \ scale);
        x = xlambda(1:n);
    else
        x = R \ (R' \ (rho*(z - u) - q));
    end          

    % z-update with relaxation
    zold = z;
    
    x_hat = alpha*x +(1-alpha)*zold;
    
    if lambdaSparse > 0
        z = shrinkage(x_hat + u, lambdaSparse/rho);
    else 
        z = min(ub, max(lb, x_hat + u));
    end
               
    % u-update
    u = u + (x_hat - z);
    
    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(P, q, r, x);

    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));
    
    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
    
    if ~QUIET
        fprintf('%3d\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end
    
    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
        break;
    end
    
    admmIters = admmIters + 1;
end    
end
