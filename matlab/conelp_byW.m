function lambda = conelp_byW(scaling,z,dims,LINSOLVER)
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
lambda(1:dims.l) = z(1:dims.l) ./ scaling.l.wl;

% Second-oder cone
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx);
    
    % multiplication
    switch( LINSOLVER )
        case 'rank1updates'
            a = scaling.q(k).a;
            zk1 = zk(2:end);
            zeta = scaling.q(k).q' * zk1;
            lambda0 = a*zk(1) - zeta;
            lambda1 = zk1 + (-zk(1) + zeta/(1+a)).*scaling.q(k).q;
            lambda(coneidx,1) = [lambda0; lambda1]./scaling.q(k).etasqrt;
            
        case 'backslash'
            lambda(coneidx,1) = scaling.q(k).W \?zk;
            
        otherwise, error('Unknown linear solver');
    end
end