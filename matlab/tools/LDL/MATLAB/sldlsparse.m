function [arg1, arg2, arg3, arg4] = ldlsparse (A, P, S, eps, delta)			    %#ok
%LDLSPARSE LDL' factorization of a real, sparse, symmetric matrix with
%regularization.
%
% Example:
%       [L, D, Parent, fl] = sldlsparse (A, P, S)
%       [L, D, Parent, fl] = sldlsparse (A, P, S, eps, delta)
%       [x, fl] = sldlsparse (A, [ ], S, b)
%       [x, fl] = sldlsparse (A, P, S, b)
%
% Let I = speye (size (A,1)). The factorization is (L+I)*D*(L+I)' = A or A(P,P).
% A must be sparse, square, and real.  Only the diagonal and upper triangular
% part of A or A(P,P) are accessed.  L is lower triangular with unit diagonal,
% but the diagonal is not returned.  D is a diagonal sparse matrix.  P is either
% a permutation of 1:n, or an empty vector, where n = size (A,1).  If not
% present, or empty, then P=1:n is assumed.  Parent is the elimination tree of
% A or A(P,P).  If positive, fl is the floating point operation count, or
% negative if any entry on the diagonal of D is zero.
%
% In the x = sldlsparse (A, P, b) usage, the LDL' factorization is not returned.
% Instead, the system A*x=b is solved for x, where both b and x are dense.
%
% REGULARIZATION:
% ---------------
% If an entry on the diagonal of D is encountered such that D[i]<=S[i]*eps,
% where S is a (permuted) sign vector, the entry is set to DELTA and the LDL'
% factorization is continued. This is called dynamic regularization, and
% this functionality is the main difference to the standard LDLSPARSE code.
%
% See also ldlsparse
%
% Copyright 2006-2007 by Timothy A. Davis, http://www.suitesparse.com
% Modified by Alexander Domahidi, Automatic Control Laboratory, ETH Zurich.

error ('sldlsparse mexFunction not found') ;
