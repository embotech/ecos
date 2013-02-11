function G = conelp_unstretch(Gtilde, dims)
% Undos the effect of conelp_stretch.
%
% Usage: G = conelp_unstretch(Gtilde, dims)
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

G = Gtilde(1:dims.l,:);
for k = 1:length(dims.q)
    coneidx = dims.l + sum(dims.q(1:k-1)) + k-1 + (1:dims.q(k));
    G = [G; Gtilde(coneidx,:)];
end
