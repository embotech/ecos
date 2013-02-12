function [L,D] = conelp_factor(K,P)
% LDL Factorization routine for CONELP solver.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% make K sparse
if( ~issparse(K) ), K = sparse(K); end

% factor & permute  
[L,D] = ldlsparse(K, P);
L = L + eye(size(L));