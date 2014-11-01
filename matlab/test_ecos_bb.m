A = [ 2., 1.
      3., 4.];

b = [4., 12.]';
c = [-1, -1]';
opts.bool_vars_idx = 2;

dims.l = 2;

[y_d, s_d, info_dd, z_d, x_d] = ecos(c, sparse(A), b, dims, opts)
