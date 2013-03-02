function [L,D,PL,QL] = conelp_lowrankfactor(A,P,S,V)
% Returns the representations of Cholesky factors from rank1 updates.
%
%   [L,D,PL,QL] = CONELP_LOWRANKFACTOR(A,P,S,V) computes L,D,PL,QL such that
%   A + V*S*V' = L*L(PL,QL)*D*L(PL,QL)'*L'. A is a matrix of size [n x n],
%   V is a [n x k] matrix with k <= n, and S is a [k x k] diagonal matrix.
%
%                    [ 1                                 ]
%                    [ p2*q1   1                         ]
%                    [ p3*q1  p3*q2    1                 ]
%          L(p,q) =  [   .      .           .            ]
%                    [   .      .               .        ]
%                    [   .      .                   .    ]
%                    [ pn*q1  pn*q2   . . . pn*qn-1   1  ]
%
%   Note that L(QL,PL) is never explicitly formed; instead, only the
%   necessary vectors p and q are returned as columns of PL and QL.
%
%   If a permutation vector P is given, it is used to determine the
%   permuted factors.
%
%   The method to compute these matrices is described in the paper by
%   D. Goldfarb and K. Scheinberg: "Product-form Cholesky factorization in
%   interior point methods for second-order cone programming", Jnl. Math.
%   Programming, 2005, p. 162. DOI: 10.1007/s10107-004-0556-1
%
% See also conelp_lowranksolve conelp_lowrankL
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

n = size(A,1);
k = size(V,2);

if( nargin < 2 )
    P = 1:n;
end

[L,D] = ldlsparse(sparse(A),P);
L = L + speye(n);

if( nargin > 2 )
    PL = NaN(n,k);
    QL = NaN(n,k);    
    s = diag(S);    
    for j = 1:k
        p = conelp_forwardsub(L,V(P,j));
        for jj = 1:j-1
            p = conelp_lowrankforwardsub(PL(:,jj),QL(:,jj),p);
        end            
        
%         Lt{j} = eye(n); % explicit representation of L(PL,QL) for debugging
        % Low-rank-update kernel, see Gill et. al., p. 516
        alpha = s(j);
        q = NaN(n,1);
        d1 = full(diag(D));
        d2 = NaN(1,n);
        for i = 1:n
            d2(i) = d1(i) + alpha*p(i)^2;
            q(i) = alpha*p(i) / d2(i);
            alpha = alpha*d1(i)/d2(i);
%             Lt{j}(i+1:end,i) = q(i)*p(i+1:end,1); % explicit representation of L(PL,QL) for debugging
        end
        
        PL(:,j) = p;
        QL(:,j) = q;
        D = sparse(diag(d2));
    end   
end
