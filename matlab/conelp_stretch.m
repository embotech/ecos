function Gtilde = conelp_stretch(G,dims)
% Adds zero rows to matrix G after each SOC cone.

n = size(G,2);

Gtilde = G(1:dims.l,:);

for k = 1:length(dims.q)
    coneidx = dims.l + sum(dims.q(1:k-1)) + (1:dims.q(k));
    Gtilde = [Gtilde; G(coneidx,:); zeros(1,n)];
end

Gtilde = sparse(Gtilde);