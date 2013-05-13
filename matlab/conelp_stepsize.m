function alpha = conelp_stepsize(lambda,ds,dz,dims,tau,dtau,kap,dkap)
% Step size computation for conelp solver.
%
%    ALPHA = CONELP_STEPSIZE(lambda,ds,dz,dims) returns the maximum step
%    size alpha such that
%
%       lambda + alpha*ds >= 0    and     lambda + alpha*dz >= 0
%
%    where the inequalities are with respect to the cones defined in dims.
%
% NOTE: The solver and the text above are heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program
%      solvers", March 20, 2010.
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

alpha = 2.0;

% LP cone
if( dims.l > 0 )
    rhok = ds(1:dims.l) ./ lambda(1:dims.l);
    sigmak = dz(1:dims.l) ./ lambda(1:dims.l);
    alpha = min([alpha, 1/max([0, -min(rhok), -min(sigmak), -dtau/tau, -dkap/kap])]);
end

% Second-order cone
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    lk = lambda(coneidx); dsk = ds(coneidx); dzk = dz(coneidx);
    
    % normalize
    lk1 = lk(2:end);
    lknorm = sqrt(lk(1)^2 - lk1'*lk1);
    lkbar = lk/lknorm;
    
    % chop vectors in pieces
    dsk0 = dsk(1); dsk1 = dsk(2:end);
    dzk0 = dzk(1); dzk1 = dzk(2:end);
    lkbar0 = lkbar(1);  lkbar1 = lkbar(2:end);
    lkbar_times_dsk = lkbar0*dsk0 - lkbar1'*dsk1;
    lkbar_times_dzk = lkbar0*dzk0 - lkbar1'*dzk1;
    rhok   = 1/lknorm*[ lkbar_times_dsk; dsk1 - (lkbar_times_dsk+dsk0)/(lkbar0+1)*lkbar1 ];
    sigmak = 1/lknorm*[ lkbar_times_dzk; dzk1 - (lkbar_times_dzk+dzk0)/(lkbar0+1)*lkbar1 ];
    
    % calculate maximum alphak
    rhok0 = rhok(1); rhok1 = rhok(2:end);
    sigmak0 = sigmak(1); sigmak1 = sigmak(2:end);
    alphak = 1/max([0, norm(rhok1,2)-rhok0, norm(sigmak1,2)-sigmak0]);
    
    % update minimum alpha
    alpha = min([alpha, alphak]);
end

if alpha > 1, alpha = 1; end
if alpha < 0, alpha = 0; end
