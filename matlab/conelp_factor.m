function [L,D,PL,QL,P] = conelp_factor(K,P,LINSOLVER,n,p,dims,scaling)
% LDL Factorization routine for CONELP solver.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.


QL = []; PL = [];
switch( LINSOLVER )
    case 'rank1updates'
        if( nargin == 7 )
            [S,Vrank1] = conelp_buildlowrankmatrices(scaling,n,p,dims);
        else
            S = []; Vrank1 = [];
        end
        [L,D,PL,QL] = conelp_lowrankfactor(K,P,S,Vrank1);
        
    case 'backslash'
        [L,D,P] = ldl(K,'vector');
        
    case 'ldlsparse'
        [L,D] = ldlsparse(sparse(K), P);
        L = L + eye(size(L));
        
    otherwise, error('Unknown linear solver');
end



