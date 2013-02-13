function x = conelp_lowrankforwardsub(P,Q,b)
% Returns the solution to L(Q,P)x = b.
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.
% Email: domahidi@control.ee.ethz.ch


[n, k] = size(P);


for j = 1:k
    
    gamma = 0;
    x = NaN(n,1);
    for i=1:n
        x(i) = b(i) - gamma*P(i,j);
        gamma = gamma + Q(i,j)*x(i);
    end
    b = x;
    
end