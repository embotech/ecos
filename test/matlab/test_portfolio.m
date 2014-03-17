D = diag(rand(5,1));
F = randn(5,2);
c = randn(5,1);

cvx_begin
    cvx_solver ecos
    variable x(5)
    minimize (c'*x + x'*D*x + sum_square(F'*x))
    sum(x) == 1
cvx_end
ecos_x = x;
ecos_optval = cvx_optval;
time_ecos = cvx_cputime;

cvx_begin
    variable x(5)
    minimize (c'*x + x'*D*x + sum_square(F'*x))
    sum(x) == 1
cvx_end
true_x = x;
true_optval = cvx_optval;
time_true = cvx_cputime;