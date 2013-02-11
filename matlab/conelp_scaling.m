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
        e = ones(dims.q(i)-1,1);
        I = eye(dims.q(i)-1);
%         Vk_pattern = [ 1,  e', 0  ;
%             e,  I , e  ;
%             0,  e', 1 ];
        Vk_pattern = [1, e'; e I];
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
W = diag(w_l);
Winv = diag( sqrt(z(1:dims.l)./s(1:dims.l)) );
S = -ones(dims.l,1);

pp = zeros(dims.l+sum(dims.q),length(dims.q));

% Second-order cone [1, ?4.2] - do this construction for all cones
for k = 1:length(dims.q)
    
    % get variables for current cone
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));
    zk = z(coneidx); sk = s(coneidx);
    
    % normalize
    zk1 = zk(2:end);
    sk1 = sk(2:end);
    zres = zk(1)^2 - zk1'*zk1;
    sres = sk(1)^2 - sk1'*sk1;
    
    % check that we're still in the cone
    if( zres <= 0 )
        %                 keyboard
        error('z not strictly in the second-order cone, exiting');
    end
    if( sres <= 0 )
        %                 keyboard
        error('s not strictly in the second-order cone, exiting');
    end
    
    % normalize
    znorm = sqrt(zres); zkbar = zk./znorm;
    snorm = sqrt(sres); skbar = sk./snorm;
    
    % Nesterov-Todd scaling point (normalized)
    gamma = sqrt((1+zkbar'*skbar)/2);
    wkbar = 1/(2*gamma)*(skbar+[zkbar(1); -zkbar(2:end)]);
    
    % Vk = stretched Wk^2 for faster factorization
    I = speye(dims.q(k)-1);
    eta = snorm/znorm;
    a = wkbar(1);
    q = wkbar(2:dims.q(k));
    v = q./sqrt(1+a);
    omega = q'*q;
    atilde = a^2 + omega;
    b = (1+a+omega/(1+a));
    c = sqrt(1+2/(1+a)+omega/((1+a)^2));
    alpha = c/b;
    qtilde = b*q;
    vtilde = c*q;
    
    Wk = sqrt(eta)*[a q'; q, I + v*v'];
    Wkinv = [a -q'; -q, I + v*v']./sqrt(eta);
    
    % get signs of D
    temp = -b^2 / atilde;
    Dk = []; Sk = [];
    for i=1:dims.q(k)-1
        Dk(i,1) = 1.0 + temp*q(i)^2;
        assert(abs(Dk(i)) > 0, 'Zero on diagonal.');
        Sk(i,1) = sign(Dk(i));
        temp = temp / Dk(i);
    end
%     rowprint(Dk,'Dk');
    Sk = [1; Sk; 2*sum(Sk<0)-1];
   
%     Vk       = eta*[atilde, qtilde', 0; qtilde, I, vtilde; 0, vtilde', -1] + diag(Sk)*1e-8;
%     Vk_noreg = eta*[atilde, qtilde', 0; qtilde, I, vtilde; 0, vtilde', -1];
    Vk =  eta*[atilde, qtilde'; qtilde, I];
    Vk_noreg =  eta*[atilde, qtilde'; qtilde, I];
    pp(coneidx(2:end),k) = vtilde;


%     [LL,DD] = ldlsparse(sparse(Vk_noreg));
    %Sk' - sign(diag(DD))'
    
    % check for numerical problems
%     assert( all(eig(Wkinv) > 0), 'Some eigenvalue of Wk^{-1} <= 0');
%     assert( all(eig(Wk) > 0), 'Some eigenvalue of Wk <= 0');
    assert( all(eig(Vk) > 0), 'Some eigenvalue of Vk <= 0');
    %assert( all(eig(Wk*Wk) > 0), 'Some eigenvalue of Wk^2 <= 0');
    
    
    % put into scaling matrix
    V = blkdiag(V, Vk);
    V_noreg = blkdiag(V_noreg, Vk_noreg);
    W = blkdiag(W, Wk);
    Winv = blkdiag(Winv, Wkinv);
%     S = [S; -Sk];
    
    % save stuff needed to apply scaling later in linear time
    scaling.q(k).q = q;
    scaling.q(k).v = v;
    scaling.q(k).qtilde = qtilde;
    scaling.q(k).vtilde = vtilde;
    scaling.q(k).a = a;
    scaling.q(k).atilde = atilde;
    scaling.q(k).alpha = alpha;
    scaling.q(k).etasqrt = sqrt(eta);
    scaling.q(k).eta = eta;
    scaling.q(k).W = Wk;
    scaling.q(k).Winv = Wkinv;
end

% calculate scaled points
lambda = conelp_timesW(scaling,z,dims); % = W*z

%% output
if( nargout >= 1 ), varargout{1} = scaling; end
if( nargout >= 2 ), varargout{2} = lambda;  end
if( nargout >= 3 ), varargout{3} = V;       end
if( nargout >= 4 ), varargout{4} = V_noreg; end
if( nargout >= 5 ), varargout{5} = pp;      end
% if( nargout >= 6 ), varargout{6} = W;       end
% if( nargout >= 7 ), varargout{7} = Winv;    end


%% DEPRECATED
%
% [p, n] = size(A);
%
% out = [scaling,lambda,V,W,S,Winv,V_noreg,D,ppp,Vpattern]
%
%     switch( LINSOLVER )
%         case 'ldlsparse5'
%             eps = 1e-3*ones(size(vtilde));
%             Vk = eta*[atilde, qtilde',    0,     0;
%                 qtilde,    I,    vtilde, eps;
%                 0,   vtilde',    -1,     0;
%                 0,     eps',      0,    -1];
%
%             o = zeros(size(eps));
%             Vk_noreg = eta*[atilde, qtilde',   0,    0;
%                 qtilde,     I,   vtilde, o;
%                 0,    vtilde',   -1,  0;
%                 0,       o',      0, -1];
%
%         case 'ldlsparse6'
%
%             Vk = eta*[atilde, qtilde';
%                 qtilde,    I];
%
%             Vk_noreg = eta*[atilde,   qtilde';
%                 qtilde,   I + vtilde*vtilde'];
%
%             newppp = [zeros(n+p+dims.l,1); zeros(sum(dims.q(1:k-1)),1); 0; sqrt(eta)*vtilde; zeros(sum(dims.q(k+1:end)),1)];
%             assert(size(newppp,1) == n+p+sum(dims.q)+dims.l);
%             ppp = [ppp, newppp];
%
%         otherwise
%
%             %         fprintf('a=%6.4e;  b=%6.4e;  c=%6.4e;  eta=%6.4e;\n',atilde,b,c,eta);
%             %         fprintf('q=[');
%             %         for ijk = 1:length(q)
%             %             fprintf('%6.4e  ',q(ijk));
%             %         end
%             %         fprintf(']'';\n');







%     end

%             % new reg method
%             temp = -eta*b^2 / atilde;
%             Dk = []; Sk = [];
%             for i=1:dims.q(k)-1
%                 Dk(i,1) = eta + temp*q(i)^2;
%                 assert(abs(Dk(i)) > 0, 'Zero on diagonal.');
%                 Sk(i,1) = sign(Dk(i));
%                 temp = temp*eta / Dk(i);
%             end
%             %             Sk = [1; Sk; 1];
%             dend = min([abs(Dk);eta*atilde]);
%             Dk = min([1e8*ones(size(Dk)), Dk],[],2);
%             Dk = sqrt(abs(blkdiag(blkdiag((eta*atilde),diag(Dk)),dend)));




