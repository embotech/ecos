function LPnetlibUFget( number,scale,wait,Method,atol )

%    LPnetlibUFget( number,scale,wait,Method,atol )
% loads LPnetlib problem "number" with UFget.
% Currently the LPnetlib problem numbers are 594--731.
%
%    LPnetlibUFget( 596,1,1 )
% loads lp_adlittle.mat, scales it, and waits inside pdco
% for possible parameter changes.
%
% The LP problem is of the form
%   min c'x  st. Ax = b, bl < x < bu.
%
%    scale  = 0  suppresses scaling (OK if A,b,c are well scaled)
%    scale  = 1  requests scaling (default)
%    wait   = 0  prevents pdco from waiting (default)
%    wait   = 1  asks pdco to wait to allow parameters to be reset.
%    Method   can be 0, 1, 2, 3, 4 or other values described in pdcoSet.m.
%    Method = 1  specifies sparse Cholesky on the K1 system AD^2A'dy = rhs.
%    Method = 3  specifies LSMR on an equivalent least-squares system.
%    Method = 4  specifies MINRES on K1.
%    Method = 21 specifies sparse LU on the SQD (symmetric quasi-definite)
%                augmented system K2.
%    atol        specifies the convergence tolerance for LSMR and MINRES
%                (default value 1e-10).

%-----------------------------------------------------------------------
% 03 Jul 2006: LPnetlib.m derived from LPnetlib3.m, which was
%              developed by Michael Saunders (MAS) and Leo Tenenblat
%              (SOL, Stanford University)
%              for experimenting with the Zoom strategy.
%              LPnetlib.m solves a single LP, with scaling as an option.
%
% 04 Jul 2006: Solution unscaled at the end.
% 04 Apr 2012: LPnetlibUFget derived from LPnetlib by Santiago Akle (SA)
%              ICME, Stanford University.
%              It uses UFget to download LPnetlib problems from
%              the University of Florida Sparse Matrix Collection
%              maintained by Tim Davis; see
%                http://www.cise.ufl.edu/research/sparse/matrices/LPnetlib/
%                http://www.cise.ufl.edu/research/sparse/mat/UFget.html
%              Currently the LPnetlib problem numbers are 594--731.
% 05 Apr 2012: (MAS) Inputs more fully documented.
% 13 Nov 2013: (MAS) Some more documentation.
%-----------------------------------------------------------------------

  if nargin < 2
     scale = 1;    % Default is to request scaling
  end
  if nargin < 3
     wait  = 1;    % Default is for pdco to wait before solving
  end
  if nargin < 4
     Method = 1;   % Default Method is Cholesky, since A is a sparse matrix
  end
  if nargin < 5
     atol = 1e-10; % Fairly accurate solve by LSQR or MINRES
  end              % if Method = 3 or 4 resp.

  % Call UFget to downlod the problem.
  Problem = UFget(number);

  A       = Problem.A;
  b       = Problem.b;
  bl      = Problem.aux.lo;
  bu      = Problem.aux.hi;
  c       = Problem.aux.c;
  [m,n]   = size(A);

  if scale
    iprint  = 1;
    scltol  = 0.9;
    [cscale,rscale] = gmscale(A,iprint,scltol);

    C = spdiags(cscale,0,n,n);   Cinv = spdiags(1./cscale,0,n,n);
    R = spdiags(rscale,0,m,m);   Rinv = spdiags(1./rscale,0,m,m);
  
    A  = Rinv*A*Cinv;
    b  = b ./rscale;
    c  = c ./cscale;
    bl = bl.*cscale;
    bu = bu.*cscale;
  end

  fixed   = find(bl==bu);
  blpos   = find(bl< bu & bl>0);
  buneg   = find(bl< bu & bu<0);
  rhs     = b - A(:,fixed)*bl(fixed) ...
              - A(:,blpos)*bl(blpos) ...
              - A(:,buneg)*bu(buneg);

  bscale  = norm(rhs,inf);   bscale  = max(bscale,1);
  oscale  = norm(c,inf);     oscale  = max(oscale,1);

  if scale
    b       = b /bscale;
    bl      = bl/bscale;
    bu      = bu/bscale;
    c       = c /oscale;
    fprintf('\n\n  Final b and c scales:  %11.4e     %11.4e', bscale, oscale)
  end
  
  
  c10     = zeros(n,1);
  c20     = zeros(n,1);
  d0      = zeros(m,1);
  gamma   = 1e-3;           % Primal regularization
  delta   = 1e-3;           % 1e-3 or 1e-4 for LP;  1 for Least squares.
  d1      = gamma;          % Can be scalar if D1 = d1*I.
  d2      = delta*ones(m,1);

%-----------------------------------------------------------------------
% Solve LP.
%-----------------------------------------------------------------------
  disp(' ');  disp(' ');  disp(' ');  disp(' ');  disp(' ');  disp(' ');
  disp('==============================================================')
  x0      = zeros(n,1);      % Initial x
  x0      = max(x0,bl);
  x0      = min(x0,bu);
  y0      = zeros(m,1);      % Initial y
  z0      = zeros(n,1);      % Initial z

  if scale
    xsize   = 1;               % Estimate of norm(x,inf) at solution
    zsize   = 1;               % Estimate of norm(z,inf) at solution
  else
    xsize   = bscale;
    zsize   = oscale;
  end

  options           = pdcoSet;
  options.mu0       = 1e-0;  % An absolute value
  options.Method    = Method;     % 1=chol  2=QR  3=LSQR  4=Minres
  options.LSMRatol1 = atol; % For LPs, LSQR must solve quite accurately
  options.LSMRMaxIter = 100; % Default = 10 (*m)
  options.OptTol    = 1e-6;
  options.FeaTol    = 1e-6;
  options.wait      = wait;  % 1 allows options to be reviewed before solve
  options.MaxIter   = 100;

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( c,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );

  if scale          % unscale x,y,z
    b  = b *bscale;
    bl = bl*bscale;
    bu = bu*bscale;
    c  = c *oscale;
    x  = x *bscale;
    y  = y *oscale;
    z  = z *oscale;

    A  = R*A*C;
    b  = b .*rscale;
    c  = c .*cscale;
    bl = bl./cscale;
    bu = bu./cscale;
    x  = x ./cscale;
    y  = y ./rscale;
    z  = z ./cscale;
    linobj = c'*x;
    fprintf('\n Unscaled linear objective value = %12.7e\n', linobj)
  else
    linobj = c'*x;
    fprintf('\n Linear objective value = %12.7e\n', linobj)
  end

  if wait
      disp('Waiting in LPnetlibUFget in case you want to look at the solution')
      keyboard
  end

