function z = conelp_kringel(x,y,dims)
% Implements the "o"-operator for conic programming.
%
% z = conelp_kringel(x,y,dims) returns z = x o y
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
z = x(1:dims.l).*y(1:dims.l);

% Second-order cone
for k = 1:length(dims.q)
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    xk = x(coneidx); yk = y(coneidx);
    clear zk
    zk(1) = xk'*yk;
    zk(2:dims.q(k),1) = xk(1).*yk(2:end) + yk(1).*xk(2:end);
    z(coneidx,1) = zk;
end
