function x = conelp_lowranksolve(L,D,P,Q,b)
% Solves the linear system 
% 
%   L*L(p1,q1)*...*L(pk,qk)*D*L(pk,qk)'*...*L(p1,q1)*L*x = b.
%
% See also conelp_lowrankfactor conelp_lowrankforwardsub conelp_lowrankbackwardsub
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.
% Email: domahidi@control.ee.ethz.ch

% forwardsubs
x = conelp_forwardsub(L,b);
x = conelp_lowrankforwardsub(P,Q,x);

% diagonal
x = conelp_byDiag(D,x);

% backwardsubs
x = conelp_lowrankbackwardsub(P,Q,x);
x = conelp_backwardsub(L',x);
