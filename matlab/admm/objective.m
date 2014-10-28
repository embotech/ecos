function obj = objective(P, q, r, x)
    obj = 0.5*x'*P*x + q'*x + r;
end
