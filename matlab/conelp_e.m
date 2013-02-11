function e = conelp_e(dims)
% Returns the conic unit vector "e" for use in CONELP solver.
%
% NOTE: The solver and the text above are heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program 
%      solvers", March 20, 2010. 
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%  
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% LP cone
e = ones(dims.l,1);

% Second-order cone
for k = 1:length(dims.q)
    e = [e; 1; zeros(dims.q(k)-1,1)]; %#ok<AGROW>
end