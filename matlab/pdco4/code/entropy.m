function [obj,grad,hess] = entropy( x )
%        [obj,grad,hess] = entropy( x )
%        computes the objective value, gradient and diagonal Hessian
%        of a separable convex function, for use with pdsco.m.
%        This is an example objective function.

%-----------------------------------------------------------------------
% 15 May 1998: Entropy function suggested by Michael Grant.
%              obj     =  sum  x(j)*log x(j)
%              grad(j) =  1 + log x(j)
%              hess(j) =  1 / x(j)
% 19 Dec 2008: hess is now a sparse matrix.
%-----------------------------------------------------------------------

   n    = length(x);
   logx = log(x);

   obj  = sum( x.*logx );
   grad = 1 + logx;
   hess = 1 ./ x;
   hess = diag(sparse(hess));  % Turn diag into a matrix for pdco2.

%-----------------------------------------------------------------------
% End function entropy
%-----------------------------------------------------------------------
