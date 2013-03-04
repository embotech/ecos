function lambda = conelp_timesW(scaling,z,dims,LINSOLVER)
% Linear time multiplication with scaling matrix W.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% assign memory
lambda = NaN(length(z),1);

% LP cone
lambda(1:dims.l,1) = scaling.l.wl .* z(1:dims.l);

% Second-oder cone
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx);
    
    % multiplication
    switch( LINSOLVER )
        
        case 'backslash'
            lambda(coneidx,1) = scaling.q(k).W*zk;
        
        case {'rank1updates', 'ldlsparse'}
            a = scaling.q(k).a;
            zk1 = zk(2:end);
            zeta = scaling.q(k).q' * zk1;
            lambda0 = a*zk(1) + zeta;
            lambda1 = zk1 + (zk(1) + zeta/(1+a)).*scaling.q(k).q;
            lambda(coneidx,1) = scaling.q(k).eta.*[lambda0; lambda1];
         
        otherwise, error('Unknown linear solver');
    end
end