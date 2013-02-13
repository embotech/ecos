function [L,D,P,Q,Lt] = conelp_lowrankfactor(A,S,V)
% Returns the representations of Cholesky factors from rank1 updates.
%
%   [L,D,P,Q] = CONELP_RANK1FACTOR(A,S,V) computes L,D,P,Q such that
%   A + V*S*V' = L*L(P,Q)*D*L(P,Q)'*L'. A is a matrix of size [n x n], 
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
%   Note that L(Q,P) is never explicitly formed; instead, only the
%   necessary vectors p and q are returned as columns of P and Q.
%
%   The method to compute these matrices is described in the paper by
%   D. Goldfarb and K. Scheinberg: "Product-form Cholesky factorization in
%   interior point methods for second-order cone programming", Jnl. Math.
%   Programming, 2005, p. 162. DOI: 10.1007/s10107-004-0556-1 
%
% See also conelp_lowranksolve conelp_lowrankL
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

[n, k] = size(V);

[L,D] = ldlsparse(sparse(A));
L = L + speye(n);

P = NaN(n,k);
Q = NaN(n,k);
d1 = full(diag(D));
d2 = NaN(1,n);
s = diag(S);
q = NaN(n,1);
for j = 1:k
    p = conelp_forwardsub(L,V(:,j));
    Lt{j} = eye(n);
    % Low-rank-update kernel, see Gill et. al., p. 516
    alpha = s(j);
    for i = 1:n
        d2(i) = d1(i) + alpha*p(i)^2;
        q(i) = alpha*p(i) / d2(i);
        alpha = alpha*d1(i)/d2(i);
        Lt{j}(i+1:end,i) = q(i)*p(i+1:end,1);
    end
    
    P(:,j) = p;
    Q(:,j) = q;
    
end
D = sparse(diag(d2));
