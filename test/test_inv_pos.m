cvx_begin
    cvx_solver ecos
    variable x
    minimize inv_pos(x)
    1 <= x <= 2
cvx_end
ecos_x = x;
ecos_optval = cvx_optval;
time_ecos = cvx_cputime;

cvx_begin
    variable x
    minimize inv_pos(x)
    1 <= x <= 2
cvx_end
true_x = x;
true_optval = cvx_optval;
time_true = cvx_cputime;