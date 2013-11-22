% Self-dual homogeneous embedding interior point method for optimization
% over linear or second-order cones. No SDP cones supported!
%
%  [x,y,info,s,z] = ECOS(c,G,h,dims,A,b,opts)
%    Solves a pair of primal and dual cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    A*x = b
%                    s >= 0
%
%        maximize    -h'*z - b'*y
%        subject to  G'*z + A'*y + c = 0
%                    z >= 0.
%
%    The inequalities are with respect to a cone C defined as the Cartesian
%    product of N + 1 cones:
% 
%        C = C_0 x C_1 x .... x C_N x C_{N+1}.
% 
%     The first cone C_0 is the nonnegative orthant of dimension ml.
%     The next N cones are second order cones of dimension mq[0], ...,
%     mq[N-1].  The second order cone of dimension m is defined as
% 
%         { (u0, u1) in R x R^{m-1} | u0 >= ||u1||_2 }.
% 
%     Input arguments:
% 
%         c is a dense matrix of size (n,1).
% 
%         dims is a struct with the dimensions of the components of C.
%         It has two fields.
%         - dims.l = ml, the dimension of the nonnegative orthant C_0.
%           (ml >= 0.)
%         - dims.q = mq = [ mq[1], mq[2], ..., mq[N] ], a row vector of N
%           integers with the dimensions of the second order cones C_1, ...,
%           C_N.  (N >= 0 and mq[k] >= 1.)
%         The default value of dims is dims.l = size(G,2) and dims.q = [].
% 
%         G is a dense or sparse matrix of size (K,n), where
% 
%             K = ml + mq[1] + ... + mq[N].
% 
%         Each column of G describes a vector
% 
%             v = ( v_0, v_1, ..., v_N+1 )
% 
%         in V = R^ml x R^mq[1] x ... x R^mq[N]
%         stored as a column vector
% 
%             [ v_0; v_1; ...; v_N+1 ].
% 
%         h is a dense matrix of size (K,1), representing a vector in V,
%         in the same format as the columns of G.
% 
%         A is a dense or sparse matrix of size (p,n).  The default value
%         is [].
% 
%         b is a dense matrix of size (p,1).   The default value is [].
%
%         opts is a struct with a single field
%         - opts.verbose, whether to be verbose (defaults to true)
% 
%         It is assumed that rank(A) = p and rank([A; G]) = n.
% 
% NOTE: The solver and the text above are heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program 
%      solvers", March 20, 2010. 
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%  
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.