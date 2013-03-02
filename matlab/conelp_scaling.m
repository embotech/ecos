function varargout = conelp_scaling(varargin)
% Computes a scaling for the CONELP solver.
%
% V = CONELP_SCALING(DIMS, LINSOLVER) returns a matrix that has the same
% sparsity pattern as a true scaling matrix, but has only ones and zeros as
% entries. This is needed to determine the ordering of the KKT matrix based
% on sparsity and symbolic factorization prior to numerical factorization.
%
% [S,LAM,V,VNR] = CONELP_SCALING(S,Z,DIMS,LINSOLVER,EPS) returns a struct S
% that contains an O(n) representation of the scaling matrix W as follows:
%
%   W = S.eta*[S.a S.q'; S.q, eye(length(S.q)) + S.b*(S.q * S.q')];
%
% or
%
%   W = D + uu' - vv',
%
% depending on the linear solver LINSOLVER chosen. LAM is the scaled
% variable LAM = W*z = W^{-1}*s, and V = W^2 + EPS*diag(S), i.e. the
% scaled scaling matrix which is regularized by a diagonal regularization
% of the right sign. VNR = W^2, i.e. the true matrix without
% regularization. This is needed in the iterative refinement steps of the
% solver.
%
% See also conelp_solve conelp_kktmatrix conelp_init conelp_timesW conelp_byW
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012-13.


%% Generate scaling pattern only
if( ~(nargin == 2 || nargin == 5) )
    error('conelp_scaling needs 2 or 5 arguments, see help conelp_scaling');
end

dims = varargin{1};
LINSOLVER = varargin{2};

if( nargin == 2 )
    
    % LP cone
    Vpattern = eye(dims.l);
    
    % SOC cone
    for i = 1:length(dims.q)
        
        switch( LINSOLVER )
            
            case 'backslash'
                Vk_pattern = ones(dims.q(i));
                
            case 'rank1updates'
                e = ones(dims.q(i),1);
                I = eye(dims.q(i));
                Vk_pattern = [I, e; e' -1];
                
            case 'ldlsparse'
                e = ones(dims.q(i),1);
                I = eye(dims.q(i));
                Vk_pattern = [I, e, e; e' -1, 0; e', 0, 1];
                
            otherwise, error('Unknown linear solver');
        end
        
        Vpattern = blkdiag(Vpattern,Vk_pattern);
    end
    
    varargout{1} = Vpattern;
    return;
end


%% Compute scalings cone-wise
if( nargin == 5 )
    % arguments
    s = varargin{1};
    z = varargin{2};
    dims = varargin{3};
    LINSOLVER = varargin{4};
    EPS = varargin{5};
    
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
        b = 1/(1+a);
        c = 1 + a + w / (1+a);
        d = 1 + 2/(1+a) + w/(1+a)^2;
        atilde = a^2 + w;
        D = eta^2*diag([atilde; ones(conesize-1,1)]);
        alpha = eta^2*d;
        beta = eta^2*c^2/d;
        
        % save stuff needed to apply scaling later in linear time
        scaling.q(k).eta = eta;
        scaling.q(k).a = a;
        scaling.q(k).b = b;
        scaling.q(k).q = q;
        scaling.q(k).u = [c/d; q];
        scaling.q(k).d = diag(D);
        scaling.q(k).alpha = alpha; 
        scaling.q(k).beta = beta;
        scaling.q(k).W = eta*[a q'; q, eye(length(q)) + b*(q*q')];
        scaling.q(k).V = scaling.q(k).W*scaling.q(k).W;
        
        % now build different scalings depending on linear solver chosen
        switch( LINSOLVER )
            
            case 'backslash'
                Vk = scaling.q(k).V + EPS*eye(conesize);
                Vk_noreg = scaling.q(k).V;
                           
            case 'rank1updates'                
                v = [1; zeros(conesize-1,1)];
                u = [c/d; q];
                Vk = [D, sqrt(alpha)*u; sqrt(alpha)*u', -1] + EPS*eye(conesize+1);
                %     Vk = D + alpha*(u*u');
                Vk_noreg = D + alpha*(u*u') - beta*(v*v');
               
                
            case 'ldlsparse'
                u0_lower = c^2/(1/w+d);
                u0_upper = min([a^2+w, c^2/d]);
                assert(u0_lower < u0_upper,'lower-upper proof does not hold');
                u0 = sqrt( u0_lower + (u0_upper - u0_lower)/2);
                u1 = c / u0;
                v0 = 0;
                v1 = sqrt(c^2/u0^2 - d);
                assert(isreal(v1),'v1 is complex');
                d1 = a^2 + w - u0^2;
                assert(d1 > 0,'d1 is negative');                
                D = blkdiag(d1,eye(conesize-1));
                u = [u0; u1*q];
                v = [v0; v1*q];
                S = [ones(conesize,1); -1; 1];
                Vk = eta^2*[D, u, v; u', -1, 0; v', 0, 1] + EPS*diag(S);
                Vk_noreg = eta^2*(D + u*u' - v*v');
                assert(all(eig(D-v*v'))>0,'D-vv'' has negative eigenvalues');
                
                
            otherwise, error('Unknown linear solver');
        end
        
        
        
        % put into scaling matrix
        V = blkdiag(V, Vk);
        V_noreg = blkdiag(V_noreg, Vk_noreg);
        
        
      
        
    end
    
    % calculate scaled points
    lambda = conelp_timesW(scaling,z,dims,LINSOLVER); % = W*z
    
    %% output
    if( nargout >= 1 ), varargout{1} = scaling; end
    if( nargout >= 2 ), varargout{2} = lambda;  end
    if( nargout >= 3 ), varargout{3} = V;       end
    if( nargout >= 4 ), varargout{4} = V_noreg; end
end

