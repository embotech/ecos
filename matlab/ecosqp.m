function [X,fval,exitflag,output,lambda,Tsolve,c,G,h,dims,Aeq,beq] = ecosqp(H,f,A,B,Aeq,Beq,lb,ub,opts)
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
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2014.

if( ~exist('opts','var') )
    opts.verbose = 1;
end

if( opts.verbose > 0 )
    disp('ECOSQP: Converting QP to SOCP...');
end

%% check if Cholesky decomposition of H exists.
assert( ~isempty(H),'Quadratic programming requires a Hessian.');
try
    W = chol(H,'upper');
catch
    warning('Hessian not positive definite, using sqrt(H) instead of chol');
    W = sqrt(H);
    k = 0;
    eliminateIdx = [];
    for i = 1:size(W,1)
        if( all(W(i,:) == 0) )
            k = k+1;
            eliminateIdx(k) = i;            
        end
    end
    W(eliminateIdx,:) = [];
    if( opts.verbose > 0 ), fprintf('%d zero rows in square root of Hessian eliminated\n',k-1); end
end
%% dimension
n = max([size(H,2), length(f)]);

%% set up SOCP problem
% The new variable is stacked as [x, a, b, t]
c = [f; 0; 0; 1];

% upper bounds
if( exist('ub','var') && ~isempty(ub) ) 
    % find indices which are upper bounded
    Aub = speye(n);
    Aub(isinf(ub),:) = [];
    Bub = ub( ~isinf(ub) );
    A = [A; Aub]; 
    B = [B; Bub];
end

% lower bounds
if( exist('lb','var') && ~isempty(lb) ) 
    % find indices which are lower bounded
    Alb = speye(n);
    Alb(isinf(lb),:) = [];
    Blb = -lb( ~isinf(lb) );
    A = [A; Alb]; 
    B = [B; Blb];
end

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
hquad = zeros(size(W,1)+2,1);
if( isempty(A) )
    G = Gquad;
    h = hquad;
    dims.l = 0;
    dims.q = size(W,1)+2;
else
    G = [A, zeros(size(A,1),3); Gquad];
    h = [B; hquad];
    dims.l = size(A,1);
    dims.q = size(W,1)+2;
end

%% sparsify
G = sparse(G);
Aeq = sparse(Aeq);

%% solve
if( opts.verbose > 0 ), fprintf('Conversion completed. Calling ECOS...\n'); end
[x,y,info,~,z] = ecos(c,G,h,dims,Aeq,beq,opts);


%% prepare return variables
X = x(1:n);
fval = c'*x;
switch( info.exitflag )
    case 1, exitflag = -2;
    case 2, exitflag = -3;
    case 0, exitflag = 1;
    otherwise, exitflag = -100;
end
output.statusstring = info.infostring;
output.iterations = info.iter;
output.time = info.timing.runtime;
lambda = [z(1:dims.l); y]; 
Tsolve = info.timing.runtime;

