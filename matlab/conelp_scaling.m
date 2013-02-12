function varargout = conelp_scaling(varargin)
% Computes a scaling for the CONELP solver.
%
% NOTE: The solver is heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program
%      solvers", March 20, 2010.
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.


%% Generate scaling pattern only
if( nargin == 1 )
    dims = varargin{1};
    
    % LP cone
    Vpattern = eye(dims.l);
    
    % SOC cone
    for i = 1:length(dims.q)
        e = ones(dims.q(i),1);
        I = eye(dims.q(i));
        Vk_pattern = [I, e; e' -1];
        Vpattern = blkdiag(Vpattern,Vk_pattern);
    end
    
    varargout{1} = Vpattern;
    return;
end

%% Compute scalings cone-wise

% arguments
s = varargin{1};
z = varargin{2};
dims = varargin{3};

% LP-cone [1, ?4.1]
w_l = sqrt( s(1:dims.l) ./ z(1:dims.l) );
V = spdiags(w_l.^2,0,dims.l,dims.l);
V_noreg = V;
scaling.l.wl = w_l;
% W = diag(w_l);
% Winv = diag( sqrt(z(1:dims.l)./s(1:dims.l)) );

% Second-order cone [1, ?4.2] - do this construction for all cones
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx); sk = s(coneidx);
    conesize = length(coneidx);
    
    % get scaling matrix and beta for rank 1 update
    s0 = sk(1); s1 = sk(2:end);
    z0 = zk(1); z1 = zk(2:end);
    
    sres = s0^2 - s1'*s1;
    zres = z0^2 - z1'*z1;
    
    assert( sres > 0, 's not in second-order cone');
    assert( zres > 0, 'z not in second-order cone');
    
    
    %% scalings
    sbar = sk./sqrt(sres);
    zbar = zk./sqrt(zres);
    eta = (sres / zres)^(1/4);
    gamma = sqrt( (1+sbar'*zbar)/2 );
    wbar = 1/2/gamma*(sbar + [zbar(1); -zbar(2:end)]);
    q = wbar(2:end);
    w = q'*q;
    a = wbar(1);
    atilde = a^2 + w;
    c = 1 + a + w / (1+a);
    d = 1 + 2/(1+a) + w/(1+a)^2;
    D = eta^2*diag([atilde; ones(conesize-1,1)]);
    alpha = eta^2*d;
    beta = eta^2*c^2/d;
    v = [1; zeros(conesize-1,1)];
    u = [c/d; q];
    Vk = [D, sqrt(alpha)*u; sqrt(alpha)*u', -1] + 1e-5*eye(length(u)+1);
%     Vk = D + alpha*(u*u');
    Vk_noreg = D + alpha*(u*u') - beta*(v*v');
    
        
    % put into scaling matrix
    V = blkdiag(V, Vk);
    V_noreg = blkdiag(V_noreg, Vk_noreg);
        
    % save stuff needed to apply scaling later in linear time
    scaling.q(k).eta = eta;
    scaling.q(k).a = a;
    scaling.q(k).q = q;
    scaling.q(k).d = diag(D);
    scaling.q(k).u = u;
    scaling.q(k).alpha = alpha;
    scaling.q(k).beta = beta;
    
end

% calculate scaled points
lambda = conelp_timesW(scaling,z,dims); % = W*z

%% output
if( nargout >= 1 ), varargout{1} = scaling; end
if( nargout >= 2 ), varargout{2} = lambda;  end
if( nargout >= 3 ), varargout{3} = V;       end
if( nargout >= 4 ), varargout{4} = V_noreg; end

