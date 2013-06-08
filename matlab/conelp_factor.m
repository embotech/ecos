function [L,D,PL,QL,P] = conelp_factor(K,P,LINSOLVER,n,p,dims,scaling, c,b,h,kap,tau)
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
        S = [ones(n,1); -ones(p,1); -ones(dims.l,1)];
        for k = 1:length(dims.q)
            S = [S; -ones(dims.q(k),1); -1; 1];
        end
        [L,D] = sldlsparse(sparse(K), P, S, 1e-14, 7e-7);
        L = L + eye(size(L));
        
%         keyboard
        
    otherwise, error('Unknown linear solver');
end



