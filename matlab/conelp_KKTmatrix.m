function K = conelp_KKTmatrix(A, G, V, eps)
% Returns regularized KKT matrix.

p = size(A,1);
[m, n] = size(G);


% eps = 0;

K = [eps*eye(n),       A',        G'     ;
          A    , -eps*eye(p), zeros(p,m) ;
          G,     zeros(m,p),     -V     ];