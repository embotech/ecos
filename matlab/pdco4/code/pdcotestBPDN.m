function pdcotestBPDN( m,n,k,lambda )

% m=50;  n=100;  k = 10;  lambda=1e-3;  pdcotestBPDN( m,n,k,lambda );
% Generates a random m by n basis pursuit denoising problem (BPDN)
% with k non zeros in the optimal solution, and treats the constraint matrix
% A as an operator.  (We need k <= n.)
%
% BPDN is the problem
%    minimize lambda ||x||_1 + 1/2 ||Ax-b||^2
% where A is m by n with m<n.
%
% The problem is tranfsormed into 
%
%    minimize    lambda e'[xl;xh] + 1/2 ||D1*[xl;xh]||^2 + 1/2 ||r||^2
%      xl,xh,r
%    subject to  [A -A]*[xl;xh] + r = b,   0 <= x <= inf,   r unconstrained,
%
% where A is m by n, e is a 2n-vector of ones,
% D1 is a small positive-definite diagonal matrix,
% and everything relevant is the same as in pdco.m.
% Here, D1 is defined by d1 in the private function toydata.

%-----------------------------------------------------------------------
% 09 Apr 2012: Example BPDN test program for pdco.m,
%              derived from pdcotestLS.m and inspired by
%              Michael Friedlander and Ewout van den Berg's spgl1 examples.
%              Santiago Akle and Michael Saunders, ICME, Stanford University.
%-----------------------------------------------------------------------

  if nargin < 4
     lambda = 1e-3;
  end

  [A,b,bl,bu,c,d1,d2,J,x_Sol] = toydata( m,n,k,lambda ); % Private function below
  
  options = pdcoSet;

  x0    = 0.5*ones(2*n,1);   % Initial x
  y0    = zeros(m,1);        % Initial y
  z0    = ones(2*n,1);       % Initial z
  xsize = 1;                 % Estimate of norm(x,inf) at solution
  zsize = 1;                 % Estimate of norm(z,inf) at solution

  options.mu0       = 1e-0;  % An absolute value
  options.LSMRatol1 = 1e-9;  % LSMR and MINRES must solve quite accurately
  options.wait      = 1;     % Allow options to be reviewed before solve

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );
  
% Rebuild the solution
  xSparse  = x(1:n) - x(n+1:end);
  figure
  hold on
  stem(xSparse,'rx');
  stem(x_Sol  ,'bo');
  hold off

% keyboard               % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% end function pdcotestBPDN
%-----------------------------------------------------------------------


function [Aop,b,bl,bu,c,d1,d2,J,xSol] = toydata( m,n,k,lambda )

%        [Aop,b,bl,bu,c,d1,d2] = toydata( m,n,k,lambda );
% defines an m by n matrix A, rhs vector b, and scalar lambda
% for the BPDN problem
%
%        minimize lambda ||x||_1 + 1/2 ||Ax-b||^2,
%
% where k is the expected number of nonzeros in x.

%-----------------------------------------------------------------------
% 10 Apr 2012: Adapted from version in pdcotestLS.m 
%-----------------------------------------------------------------------

  randn('state',0);
  [Q,R]   = qr(randn(n,m),0); % Form a dense matrix with orthogonal rows
  A       = Q';

  p       = randperm(n);      % Pick a random support  
  xSol    = zeros(n,1);       % Form a sparse solution
  J       = p(1:k);
  xSol(J) = randn(k,1);
  b       = A*xSol;           % Form the corresponding rhs

  % Put the data into the format that PDCO understands
  Aop     = @(mode,m,n,x) operatorA(mode,x,A);
  c       =  ones(2*n,1)*lambda;
  bl      = zeros(2*n,1);     % x > 0
  bu      =   inf(2*n,1);
  d1      = 1e-10;            % Primal regularization
  d2      = 1;                % Make it a regularized least squares problem
%-----------------------------------------------------------------------
% end private function toydata
%-----------------------------------------------------------------------


function y = operatorA( mode,x,A )

% Computes products with [A -A] and [A -A]' to be used
% by pdco as an operator.

   [m,n] = size(A);
   if mode == 1
      y = A*(x(1:n) - x(n+1:2*n));
   else
      r = A'*x;
      y = [r;-r];
   end

%----------------------------------------------------------------------
% End private function operatorA
%----------------------------------------------------------------------
