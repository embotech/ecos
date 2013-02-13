function x = conelp_lowrankbackwardsub(P,Q,b)
% Returns the solution to L(Q,P)^T x = b.

[n, k] = size(P);
x = NaN(n,1);

for j = 1:k
    
    gamma = 0;
    for i=n:-1:1
        x(i) = b(i) - gamma*Q(i,j);
        gamma = gamma + P(i,j)*x(i);
    end
    b = x;
    
end
