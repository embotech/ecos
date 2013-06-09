function [P,xhat,yhat,shat,zhat,kaphat,tauhat] = conelp_init(c,G,Gtilde,h,dims,A,b, LINSOLVER, EPS, NITREF, LINSYSACC)
% Initialization of variables for conelp solver, see [1, ?7.3].
%
% NOTE: The solver and the text above are heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program
%      solvers", March 20, 2010.
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

%% dimensions
[mtilde, n] = size(Gtilde);
m = dims.l + sum(dims.q);
p = size(A,1);


%% compute permutation
Vpattern = conelp_scaling(dims, LINSOLVER, 'pattern');
Kpattern = conelp_KKTmatrix(A,Gtilde,Vpattern,0);
P = conelp_getPerm(Kpattern~=0);


%% assemble and factor coefficient matrix
Vinit = conelp_scaling(dims, LINSOLVER, 'init');
Kinit = conelp_KKTmatrix(A,Gtilde,Vinit,EPS);
% [Linit,Dinit] = conelp_factor(Kinit,P);
[Linit,Dinit,~,~,P] = conelp_factor(Kinit,P,LINSOLVER,n,p,dims);


%% primal variables
%  * solve xhat = arg min ||Gx-h||_2^2  such that Ax = b
%  * r = h - G*xhat
% These two equations are solved by
%
% [ 0   A'  G' ] [ xhat ]     [ 0 ]
% [ A   0   0  ] [  y   ]  =  [ b ]
% [ G   0  -I  ] [ -r   ]     [ h ]
%
% and then take shat = r if alphap < 0, zbar + (1+alphap)*e otherwise
% where alphap = inf{ alpha | sbar + alpha*e >= 0 }

% v = M\[zeros(n,1); b; h];
% xhat = v(1:n);
% r = -v(n+p+1:end);

[xhat, ~, minus_r] = conelp_solve(Linit,Dinit,P,[],[], zeros(n,1),b,h, A,G,Vinit, dims, NITREF,LINSOLVER, LINSYSACC);
shat = bring2cone(-minus_r,dims);


%% dual variables
% solve (yhat,zbar) = arg min ||z||_2^2 such that G'*z + A'*y + c = 0
%
% we can solve this by
%
% [ 0   A'  G' ] [  x   ]     [ -c ]
% [ A   0   0  ] [ yhat ]  =  [  0 ]
% [ G   0  -I  ] [ zbar ]     [  0 ]
%
% and then take zhat = zbar if alphad < 0, zbar + (1+alphad)*e otherwise
% where alphad = inf{ alpha | zbar + alpha*e >= 0 }

%v = M\[-c; zeros(p,1); zeros(m,1)];
%yhat = v(n+1:n+p);
%zbar = v(n+p+1:end);


[~, yhat, zbar] = conelp_solve(Linit,Dinit,P,[],[], -c,zeros(p,1),zeros(m,1), A,G,Vinit, dims, NITREF,LINSOLVER,LINSYSACC);
zhat = bring2cone(zbar,dims);

%% homogeneous embedding variables
kaphat = 1.0;
tauhat = 1.0;

end


%% brings slack variables or multipliers to cone
function [shat, alpha] = bring2cone(r,dims)
% checks whether r > 0 where > is w.r.t. the cone defined by dims
% alpha is the biggest residual
alpha = -1;

% LP cone
if ( any(r(1:dims.l) <= 0 ) )
    alpha = -min(r(1:dims.l));
end

% Second-order cone
for k = 1:length(dims.q)
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    rk = r(coneidx);
    res = rk(1) - norm(rk(2:end),2);
    if( res <= 0 )
        alpha = max([alpha,-res]);
    end
end

shat = r + (1+alpha)*conelp_e(dims);
end