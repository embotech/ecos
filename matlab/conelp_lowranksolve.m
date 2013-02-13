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
if( ~isempty(P) && ~isempty(Q) )
    x = conelp_lowrankforwardsub(P,Q,x);
end

% diagonal
x = conelp_byDiag(D,x);

% backwardsubs
if( ~isempty(P) && ~isempty(Q) )
    x = conelp_lowrankbackwardsub(P,Q,x);
end
x = conelp_backwardsub(L',x);
