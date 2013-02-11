function lambda = conelp_timesWsquare(scaling,z,dims)
% Linear time multiplication with scaling matrix W.
%
% NOTE: The solver is heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program
%      solvers", March 20, 2010.
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% assign memory
lambda = NaN(length(z),1);

% LP cone
lambda(1:dims.l) = scaling.l.wl .* scaling.l.wl .* z(1:dims.l);

% Second-oder cone
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx);
    
    % multiplication
    atilde = scaling.q(k).atilde;  
    beta = scaling.q(k).alpha^2;
    zk1 = zk(2:end);
    zeta = scaling.q(k).qtilde' * zk1;
    lambda0 = atilde*zk(1) + zeta;    
    lambda1 = zk1 + (zk(1) + zeta*beta).*scaling.q(k).qtilde;
    lambda(coneidx) = scaling.q(k).etasqrt.*[lambda0; lambda1].*scaling.q(k).etasqrt;
end