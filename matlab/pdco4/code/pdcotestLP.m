function pdcotestLP( m,n )

% m=50;  n=100;  pdcotestLP( m,n );
% Generates a random m by n LP problem
%   min c'x  st. Ax = b, bl < x < bu,
% and runs it on pdco.m.

%-----------------------------------------------------------------------
% 08 Oct 2002: Simple test program for pdco.m.
%              "A" is an explicit sparse matrix, not a function.
%              Michael Saunders, SOL, Stanford University.
% 16 Sep 2003: "c" is now an explicit vector.
%              It may be passed to pdco instead of function 'linobj'.
%-----------------------------------------------------------------------

  [A,b,bl,bu,c,d1,d2] = toydata( m,n ); % Private function below
% D  = sum(A,1);   D(find(D==0)) = 1;
% D  = sparse( 1:n, 1:n, 1./D, n, n );
% A  = A*D;                             % Normalize cols of A

  options = pdcoSet;

  x0    = ones(n,1)/n;       % Initial x
  y0    = zeros(m,1);        % Initial y
  z0    = ones(n,1)/n;       % Initial z
  xsize = 1;                 % Estimate of norm(x,inf) at solution
  zsize = 1;                 % Estimate of norm(z,inf) at solution

  options.mu0       = 1e-0;  % An absolute value
  options.LSQRatol1 = 1e-8;  % For LPs, LSQR must solve quite accurately
  options.wait      = 1;     % Allow options to be reviewed before solve

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );

% keyboard                   % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% End function pdcotestLP
%-----------------------------------------------------------------------


function [A,b,bl,bu,c,d1,d2] = toydata( m,n )

%        [A,b,bl,bu,c,d1,d2] = toydata( m,n );
%        defines an m by n matrix A, rhs vector b, and cost vector c
%        for use with pdco.m.

%-----------------------------------------------------------------------
% 12 Feb 2001: First version of toydata.m.
% 30 Sep 2002: pdsco version modified for pdco.
% 16 Sep 2003: c now generated.  It may be passed to pdco.m.
% 23 Sep 2003: Simplified a bit to match pdcotestLS.m.
%-----------------------------------------------------------------------

  rand('state',0);
  density = 0.25;
  rc      = 1e-1;

  em      = ones(m,1);
  en      = ones(n,1);
  zn      = zeros(n,1);
  bigbnd  = 1e+30;

  A       = sprand(m,n,density,rc);
  x       = en;

  gamma   = 1e-3;       % Primal regularization
  delta   = 1e-3;       % 1e-3 or 1e-4 for LP;  1 for Least squares.

  d1      = gamma;      % Can be scalar if D1 = d1*I.
  d2      = delta*em;

  b       = full(A*x);
  c       = rand(n,1);

  bl      = zn;         % x > 0
  bu      = 10*en;      % x < 10
% bl(1)   = - bigbnd;   % Test "free variable" (no bounds)
% bu(1)   =   bigbnd;

%-----------------------------------------------------------------------
% End private function toydata
%-----------------------------------------------------------------------
