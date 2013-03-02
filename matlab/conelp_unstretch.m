function G = conelp_unstretch(Gtilde, dims, N)
% Undos the effect of conelp_stretch.
%
% Usage: G = conelp_unstretch(Gtilde, dims, N)
%
% See also conelp_unstretch
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

if( N == 0 )
    G = Gtilde;
else
    
    G = Gtilde(1:dims.l,:);
    for k = 1:length(dims.q)
        coneidx = dims.l + sum(dims.q(1:k-1)) + N*(k-1) + (1:dims.q(k));
        G = [G; Gtilde(coneidx,:)];
    end
    
end
