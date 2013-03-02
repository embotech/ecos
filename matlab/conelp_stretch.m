function Gtilde = conelp_stretch(G,dims,N)
% Expand conic inequality coefficient matrix G by adding zero rows.
%
% GTILDE = CONELP_STRETCH(G,DIMS,N) adds N zero rows to matrix G after
% each SOC cone.
%
% See also conelp_unstretch
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012-13.

if( N == 0 )
    Gtilde = G;
else   
    
    n = size(G,2);    
    Gtilde = G(1:dims.l,:);
    
    for k = 1:length(dims.q)
        coneidx = dims.l + sum(dims.q(1:k-1)) + (1:dims.q(k));
        Gtilde = [Gtilde; G(coneidx,:); zeros(N,n)];
    end
    
    Gtilde = sparse(Gtilde);
end