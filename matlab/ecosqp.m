function [X,fval,exitflag,output,lambda,Tsolve] = ecosqp(H,f,A,B,Aeq,Beq,lb,ub)
%QUADPROG Quadratic programming. 
%   X = ECOSQP(H,f,A,b) attempts to solve the quadratic programming 
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b 
%             x    
%
%   X = ECOSQP(H,f,A,b,Aeq,beq) solves the problem above while 
%   additionally satisfying the equality constraints Aeq*x = beq.
%
%   X = ECOSQP(H,f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution is in the 
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds 
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if 
%   X(i) is unbounded above.
%
% This file is an interface for ECOS, and rewrites the QP as a second-order
% cone program (SOCP) that can be solved by ECOS. The rewriting occurs
% through introduction of 3 variables t,a,b and rewriting the problem as
%
%     min     f'x + t
%   x,t,a,b
%
%     s.t.    A*x <= b
%           [ lb <= x <= ub ]
%           [ Aeq*x == beq    ]
%
%              a = (t-1)/2      }
%              b = (t+1)/2      } these constraints say 0.5*x'*H*x <= t
%             ||W*x ||         } 
%             || a  ||2  <= b  }
%
% where L is such that W'*W = 0.5*H (--> L = chol(0.5*H) )
%
% See also ECOS
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

%% check if Cholesky decomposition of H exists.
assert( ~isempty(H),'Quadratic programming requires a Hessian.');
W = chol(H,'upper');

%% dimension
n = max([size(H,1), length(f)]);

%% set up SOCP problem
% The new variable is stacked as [x, a, b, t]
c = [f; 0; 0; 1];

% upper/lower bounds
if( exist('ub','var') ), A = [A; speye(n)]; B = [B; ub]; end
if( exist('lb','var') ), A = [A; -speye(n)]; B = [B; -lb]; end

% a = 0.5*t - 0.5, 
% b = 0.5*t + 0.5
Aeq_a = [zeros(1,n), 1, 0, -0.5]; beq_a = -0.5;
Aeq_b = [zeros(1,n), 0, 1, -0.5]; beq_b = 0.5;
if( isempty(Aeq) )
    Aeq = [Aeq_a; Aeq_b]; beq = [beq_a; beq_b];
else 
    Aeq = [Aeq, zeros(size(Aeq,1),3)];
    Aeq = [Aeq; Aeq_a; Aeq_b]; beq = [Beq; beq_a; beq_b];
end

% add second-order cone constraint
Gquad = -[zeros(1,n), 0, sqrt(2), 0;
             W,     zeros(size(W,1),3);
         zeros(1,n), sqrt(2), 0, 0];
hquad = zeros(n+2,1);
if( isempty(A) )
    G = Gquad;
    h = hquad;
    dims.l = 0;
    dims.q = n+2;
else
    G = [A, zeros(size(A,1),3); Gquad];
    h = [B; hquad];
    dims.l = size(A,1);
    dims.q = n+2;
end

%% sparsify
G = sparse(G);
Aeq = sparse(Aeq);

%% solve
[x,y,info,~,z] = ecos(c,G,h,dims,Aeq,beq);

%% prepare return variables
X = x(1:n);
fval = c'*x;
switch( info.exitflag )
    case 1, exitflag = -2;
    case 2, exitflag = -3;
    case 0, exitflag = 1;
end
output = info.infostring;
lambda = [z(1:dims.l); y]; 
Tsolve = info.timing.runtime;

