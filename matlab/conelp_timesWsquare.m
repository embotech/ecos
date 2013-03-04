function lambda = conelp_timesWsquare(scaling,z,dims, LINSOLVER)
% Linear time multiplication with square of scaling matrix W.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% assign memory
lambda = NaN(length(z),1);

% LP cone
lambda(1:dims.l) = scaling.l.wl .* z(1:dims.l) .* scaling.l.wl;

% Second-oder cone
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx);
    
    % multiplication
    switch( LINSOLVER )
        
        case 'backslash'
            lambda(coneidx,1) = scaling.q(k).V*zk;
        
        case 'rank1updates'
            d = scaling.q(k).d;
            u = scaling.q(k).u;
            alpha = scaling.q(k).alpha;
            beta = scaling.q(k).beta;
            lambda(coneidx,1) = d.*zk + alpha*u*(u'*zk) - [beta*zk(1); zeros(dims.q(k)-1,1)];
            
        case 'ldlsparse'
            u = scaling.q(k).u;
            v = scaling.q(k).v;
            d = scaling.q(k).d;
            uz = u'*zk;
            vz = v'*zk;
            lambda(coneidx,1) = scaling.q(k).eta^2.*(d.*zk + u.*uz - v.*vz);
            
        otherwise, error('Unknown linear solver');
            
    end
end
