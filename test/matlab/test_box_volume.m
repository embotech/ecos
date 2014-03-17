% Box volume maximization
% Boyd, Kim, Vandenberghe, and Hassibi, "A Tutorial on Geometric Programming"
% Written for CVX by Almir Mutapcic 02/08/06
% (a figure is generated)
%
% Maximizes volume of a box-shaped structure which has constraints
% on its total wall area, its total floor area, and which has lower
% and upper bounds on the aspect ratios. This leads to a GP:
%
%   maximize   h*w*d
%       s.t.   2(h*w + h*d) <= Awall, w*d <= Afloor
%              alpha <= h/w <= beta
%              gamma <= d/w <= delta
%
% where variables are the box height h, width w, and depth d.

% problem constants
alpha = 0.5; beta = 2; gamma = 0.5; delta = 2;

N = 10;
Afloor = 100;
Awall  = 1000;


% solve a GP
cvx_begin gp
  cvx_solver ecos
  variables h w d
  % objective function is the box volume
  maximize( h*w*d )
  subject to
    2*(h*w + h*d) <= Awall;
    w*d <= Afloor;
    alpha <= h/w <= beta;
    gamma <= d/w <= delta;
cvx_end
ecos_h = h;
ecos_w = w;
ecos_d = d;
ecos_optval = cvx_optval;
time_ecos = cvx_cputime;

cvx_begin gp
  variables h w d
  % objective function is the box volume
  maximize( h*w*d )
  subject to
    2*(h*w + h*d) <= Awall;
    w*d <= Afloor;
    alpha <= h/w <= beta;
    gamma <= d/w <= delta;
cvx_end
true_h = h;
true_w = w;
true_d = d;
true_optval = cvx_optval;
time_true = cvx_cputime;