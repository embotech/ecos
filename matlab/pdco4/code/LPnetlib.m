function LPnetlib(name,scale,wait)

% LPnetlib( 'lp_adlittle' );
% Loads an LP problem from lp_adlittle.mat and runs it on pdco.m.
% The LP problem is of the form
%   min c'x  st. Ax = b, bl < x < bu.
%
% LPnetlib( 'lp_adlittle',scale,wait );
%    scale = 0  suppresses scaling (OK if A,b,c are well scaled)
%    scale = 1  requests scaling (default)
%    wait  = 0  prevents pdco from waiting (default)
%    wait  = 1  asks pdco to wait to allow parameters to be reset.
%
% The *.mat file format is used by Tim Davis for the
% LPnetlib group of sparse matrices:
%    http://www.cise.ufl.edu/research/sparse/matrices/LPnetlib/
%
% Files such as 'lp_adlittle.mat' should be accessible via the
% current path.  For example,
%    addpath ~/matlab/LPnetlib
% if that directory contains some of the LPnetlib datasets.
    
%-----------------------------------------------------------------------
% 03 Jul 2006: LPnetlib.m derived from LPnetlib3.m, which was
%              developed by Michael Saunders and Leo Tenenblat
%              (SOL, Stanford University)
%              for experimenting with the Zoom strategy.
%              LPnetlib.m solves a single LP, with scaling as an option.
%
% 04 Jul 2006: Solution unscaled at the end.
% 13 Nov 2013: (San Kim) The Problem structure now matches the
%              revised format used by Tim Davis in the LPnetlib group.
%-----------------------------------------------------------------------

  if nargin < 2
     scale = 1;   % Default is to request scaling
  end
  if nargin < 3
     wait  = 1;
  end

  load(name);

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
  zoom    = 1e-3;

%-----------------------------------------------------------------------
% Solve LP accurately for comparison.
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
  options.Method    = 1;     % 1=chol  2=QR  3=LSQR
  options.LSQRatol1 = 1e-10; % For LPs, LSQR must solve quite accurately
  options.LSQRMaxIter = 100; % Default = 10 (*m)
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
    fprintf('\n Unscaled linear objective = %12.7e\n', linobj)
  end

  if wait
    disp('Waiting in LPnetlib in case you want to look at the solution')
    keyboard
  end
