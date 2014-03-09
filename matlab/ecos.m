% ECOS - Embedded COnic Solver.
%
% Self-dual homogeneous embedding interior point implementation for optimization
% over linear or second-order cones. ECOS does not support semi-definite
% cones.
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims) Solves a pair of primal and dual
%   cone programs
% 
%        minimize    c'*x
%        subject to  G*x + s = h
%                    s >= 0
%
%        maximize    -h'*z
%        subject to  G'*z + c = 0
%                    z >= 0.
%
%   The inequalities are with respect to a cone C defined as the Cartesian
%   product of N + 1 cones:
% 
%        C = C_0 x C_1 x .... x C_N x C_{N+1}.
% 
%     The first cone C_0 is the nonnegative orthant of dimension dims.l.
%     The next N cones are second order cones of dimension dims.q(1), ...,
%     dims.q(N), where a second order cone of dimension m is defined as
% 
%         { (u0, u1) in R x R^{m-1} | u0 >= ||u1||_2 }.
% 
%     Input arguments:
% 
%         c is a dense matrix of size (n,1) (column vector)
% 
%         dims is a struct with the dimensions of the components of C.
%         It has two fields.
%         - dims.l, the dimension of the nonnegative orthant C_0, with l>=0.
%         - dims.q, a row vector of N integers with the dimensions of the 
%           second order cones C_1, ..., C_N. (N >= 0 and q(1) >= 1.)
% 
%         G is a sparse matrix of size (K,n), where
% 
%             K = dims.l + dims.q(1) + ... + dims.q(N).
%               = dims.l + sum(dims.q)
% 
%         Each column of G describes a vector
% 
%             v = ( v_0, v_1, ..., v_N+1 )
% 
%         in V = R^dims.l x R^dims.q(1) x ... x R^dims.q(N)
%         stored as a column vector
% 
%             [ v_0; v_1; ...; v_N+1 ].
% 
%         h is a dense matrix of size (K,1), representing a vector in V,
%         in the same format as the columns of G.
%
%
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b) Solves a pair of primal and 
%   dual cone programs
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
%      where c,G,h,dims are defined as above, and A is a sparse matrix of 
%      size (p,n), and b is a dense matrix of size (p,1).
%       
%      It is assumed that rank(A) = p and rank([A; G]) = n.
%
% 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,otps) and 
%   [x,y,info,s,z] = ECOS(c,G,h,dims,A,b,otps) are as above, with the struct 
%   opts used to control settings of the solver. The following fields can
%   be present:
%
%      .verbose - whether to inform on progess (default: true)
%      .feastol - stopping tolerance on infeasibilities (default: 1e-5)
%      .abstol  - stopping tolerance (default: 1e-6)
%      .reltol  - rel. stopping toletance (default: 1e-6)
%      .maxit   - maximum number of iterations (default: 30)
% 
%
% Details on ECOS can be found at github.com/ifa-ethz/ecos and in the paper 
%    Alexander Domahidi, Eric Chu, Stephen Boyd. "ECOS: An Embedded
%    Conic Solver." In proceedings of European Control Conference (ECC), 
%    pp. 3071-3076, Zurich, Switzerland, July 2013."
%
% COPYING: ECOS is under GPLv3. For commercial licenses and support email 
% to ecos@embotech.com.
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2012-2014.
%
% See also ECOS_LICENSE