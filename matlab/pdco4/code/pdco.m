function [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco(pdObj,pdMat,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize)

%-----------------------------------------------------------------------
% pdco.m: Primal-Dual Barrier Method for Convex Objectives (23 Nov 2013)
%-----------------------------------------------------------------------
%        [x,y,z,inform,PDitns,CGitns,time] = ...
%   pdco(pdObj,pdMat,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize);
%
% solves optimization problems of the form
%
%    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
%
% where
%    phi(x) is a smooth convex function  defined by function pdObj;
%    A      is an m x n matrix defined by matrix or function pdMat;
%    b      is a given m-vector;
%    D1, D2 are positive-definite diagonal matrices defined from d1, d2.
%           In particular, d2 indicates the accuracy required for
%           satisfying each row of Ax = b.
%
% D1 and D2 (via d1 and d2) provide primal and dual regularization
% respectively.  They ensure that the primal and dual solutions
% (x,r) and (y,z) are unique and bounded.
%
% A scalar d1 is equivalent to d1 = ones(n,1), D1 = diag(d1).
% A scalar d2 is equivalent to d2 = ones(m,1), D2 = diag(d2).
% Typically, d1 = d2 = 1e-4.
% These values perturb phi(x) only slightly  (by about 1e-8) and request
% that A*x = b be satisfied quite accurately (to about 1e-8).
% Set d1 = 1e-4, d2 = 1 for least-squares problems with bound constraints.
% The problem is then equivalent to
%
%    minimize    phi(x) + 1/2 norm(d1*x)^2 + 1/2 norm(A*x - b)^2
%    subject to  bl <= x <= bu.
%
% More generally, d1 and d2 may be n and m vectors containing any positive
% values (preferably not too small, and typically no larger than 1).
% Bigger elements of d1 and d2 improve the stability of the solver.
%
% At an optimal solution, if x(j) is on its lower or upper bound,
% the corresponding z(j) is positive or negative respectively.
% If x(j) is between its bounds, z(j) = 0.
% If bl(j) = bu(j), x(j) is fixed at that value and z(j) may have
% either sign.
%
% Also, r and y satisfy r = D2 y, so that Ax + D2^2 y = b.
% Thus if d2(i) = 1e-4, the i-th row of Ax = b will be satisfied to
% approximately 1e-8.  This determines how large d2(i) can safely be.
%
%
% EXTERNAL FUNCTIONS:
% options         = pdcoSet;                  provided with pdco.m
% [obj,grad,hess] = pdObj( x       );         provided by user
%               y = pdMat( mode,m,n,x );      provided by user if pdMat
%                                             is a handle, not a matrix
% lsmr.m    is required if options.Method = 3.
% minres.m  is required if options.Method = 4.
%
% INPUT ARGUMENTS:
% pdObj      may be an explicit n x 1 column vector c,
%            or a function handle pdObj
%            such that [obj,grad,hess] = pdObj(x) defines
%            obj  = phi(x)              : a scalar,
%            grad = gradient of phi(x)  : an n-vector,
%            hess = diag(Hessian of phi): an n-vector.
%         Examples:
%            If phi(x) is the linear function c'*x, pdObj could be
%            be the vector c, or the handle of a function that returns
%               [obj,grad,hess] = [c'*x, c, zeros(n,1)].
%            If phi(x) is the entropy function E(x) = sum x(j) log x(j),
%            pdObj should return
%               [obj,grad,hess] = [E(x), log(x)+1, 1./x].
% pdMat      may be an explicit m x n matrix A (preferably sparse!),
%            or a function handle such that  y = pdMat(mode,m,n,x )
%            returns   y = A*x (mode=1)  or  y = A'*x (mode=2).
% b          is an m-vector.
% bl         is an n-vector of lower bounds.  Non-existent bounds
%            may be represented by bl(j) = -Inf or bl(j) <= -1e+20.
% bu         is an n-vector of upper bounds.  Non-existent bounds
%            may be represented by bu(j) =  Inf or bu(j) >=  1e+20.
% d1, d2     may be positive scalars or positive vectors (see above).
% options    is a structure that may be set and altered by pdcoSet
%            (type help pdcoSet).
% x0, y0, z0 provide an initial solution.
% xsize, zsize are estimates of the biggest x and z at the solution.
%            They are used to scale (x,y,z).  Good estimates
%            should improve the performance of the barrier method.
%
%
% OUTPUT ARGUMENTS:
% x          is the primal solution.
% y          is the dual solution associated with Ax + D2 r = b.
% z          is the dual solution associated with bl <= x <= bu.
% inform = 0 if a solution is found;
%        = 1 if too many iterations were required;
%        = 2 if the linesearch failed too often;
%        = 3 if the step lengths became too small;
%        = 4 if Cholesky said ADDA was not positive definite.
% PDitns     is the number of Primal-Dual Barrier iterations required.
% CGitns     is the number of Conjugate-Gradient  iterations required
%            if an iterative solver is used (LSMR or MINRES).
% time       is the cpu time used (via cputime).
%            We also use tic/toc to allow for multicore systems.
%----------------------------------------------------------------------

% PRIVATE FUNCTIONS:
%    pdxxxbounds
%    pdxxxdistrib
%    pdxxxlsmrmat
%    pdxxxminresmat
%    pdxxxmat
%    pdxxxmerit
%    pdxxxresid1
%    pdxxxresid2
%    pdxxxstep
%
% NOTES:
% The matrix A should be reasonably well scaled: norm(A,inf) =~ 1.
% The vector b and objective phi(x) may be of any size, but ensure that
% xsize and zsize are reasonably close to norm(x,inf) and norm(z,inf)
% at the solution.
%
%
% AUTHOR:
%    Michael Saunders, Systems Optimization Laboratory (SOL),
%    Dept of Management Science and Engineering (MS&E), and
%    Institute of Computational and Mathematical Engineering (ICME),
%    Stanford University, Stanford, CA 94305, USA (saunders@stanford.edu).
%
%    http://www.stanford.edu/group/SOL/   (SOL)
%    http://www.stanford.edu/dept/MSandE/ (MS&E)
%    http://icme.stanford.edu/            (ICME)
%    http://www.stanford.edu/~saunders/
%
% CONTRIBUTORS:
%    (BK ) Byunggyoo Kim, Samsung, Seoul, Korea     (hightree@samsung.co.kr)
%    (CMM) Chris Maes, ICME, Stanford University    (cmaes@stanford.edu)
%    (SA ) Santiago Akle, ICME, Stanford University (akle@stanford.edu)
%    (MZ ) Matt Zahr, ICME, Stanford University     (mzahr@stanford.edu)
%
% DEVELOPMENT:
% 20 Jun 1997: Original version of pdsco.m derived from pdlp0.m.
% 29 Sep 2002: Original version of pdco.m  derived from pdsco.m.
%              Introduced D1, D2 in place of gamma*I, delta*I
%              and allowed for general bounds bl <= x <= bu.
% 06 Oct 2002: Allowed for fixed variabes: bl(j) = bu(j) for any j.
% 15 Oct 2002: Eliminated some work vectors (since m, n might be LARGE).
%              Modularized residuals, linesearch
% 16 Oct 2002: pdxxx..., pdDDD... names rationalized.
%              pdAAA eliminated (global copy of A).
%              pdMat is now used directly as an explicit A or a function.
%              NOTE: If pdMat is a function, it now has an extra parameter.
% 01 Nov 2002: Bug fixed in feval in pdxxxmat.
% 20 Jan 2003: Revised for Tomlab by Kenneth Holmstrom
% 19 Apr 2003: Bug fixed in pdxxxbounds.
% 07 Aug 2003: Let d1, d2 be scalars if input that way.
% 10 Aug 2003: z isn't needed except at the end for output.
% 10 Aug 2003: mu0 is now an absolute value -- the initial mu.
% 13 Aug 2003: Access only z1(low) and z2(upp) everywhere.
%              stepxL, stepxU introduced to keep x within bounds.
%              (With poor starting points, dx may take x outside,
%              where phi(x) may not be defined.
%              Entropy once gave complex values for the gradient!)
% 16 Sep 2003: pdObj can now be a vector c, implying a linear obj c'*x.
% 19 Sep 2003: Large system K4 dv = rhs implemented.
% 23 Sep 2003: Options LSproblem and LSmethod replaced by Method.
% 18 Nov 2003: stepxL, stepxU gave trouble on lptest (see 13 Aug 2003).
%              Disabled them for now.  Nonlinear problems need good x0.
% 19 Nov 2003: Bugs with x(fix) and z(fix).
%              In particular, x(fix) = bl(fix) throughout, so Objective
%              in iteration log is correct for LPs with explicit c vector.
% 21 Nov 2003: Second Revision for Tomlab by Kenneth Holmstrom
%              Do all output if PriLev > 0, PriLev added to pdxxxbounds
% 28 Nov 2003: pdxxxmat: y = pdMat'*x; replaced by y = (x'*pdMat)';
% 12 Jan 2006: If options.mu0 = 0, initialize mu0 from the average
%              complementarity.

% 23 Apr 2008: (CMM) Converted strings pdObj and pdMat to function handles.
%              TomLab changes deleted for now.
%              Removed global variables pdDDD1, pdDDD2, pdDDD3.
%              Removed internal pdxxxlsqr routine (we call new lsqr directly).
%              Created new private function pdxxxstop.
%              Modified pdxxxmat to work with function handles. 
% 17 May 2008: Don't let mu get less than 0.1*max([Pinf Dinf Cinf]).
% 29 Aug 2008: Initialize x(pos) = x1(pos). Trying to keep x>bl (usually 0).
% 30 Mar 2009: If Method<=4, the Hessian must be diagonal as for the
%              previous PDCO, so we can apply LSQR to a least-squares
%              subproblem.
% 26 Feb 2010: Treat zero bounds specially using zlo and zup instead of pos.
%              zlo and zup now point to nonfixed variables for which
%              bl = 0 or bu = 0  (with bl < bu).
%              We keep x(zlo) = x1(zlo) and x(zup) = - x2(zup), so the
%              corresponding x is always strictly positive or negative.
% 03 Apr 2010: Ronan Fleming noticed the output z = grad - A'y.
%              It now includes (d1.^2).*x from the primal regularization.
% 03 Apr 2010: If all variables are free, define center = 1.0.
% 02 Jun 2011: options.backtrack added.  The default value 0
%              turns the backtracking linesearch OFF
%              because it seems to inhibit the convergence of some problems
%              arising in a fixed-point iteration for systems biology.
%              options.backtrack = 1 will backtrack if necessary to
%              reduce the 2-norm of the residuals of the Newton system
%              for computing a search direction (dx,dy,dz).
% 03 Apr 2012: (SA) PDCO now uses LSMR when Method=3 (in place of LSQR).
%              Removed pdxxxstop, renamed pdxxxlsqrmat to pdxxxlsmrmat.
%              PDCO now uses MINRES when Method=4.
%              pdxxxminresmat defines the MINRES operator.
% 28 Apr 2012: Added diag preconditioning for MINRES when pdMat is a matrix
%              (similar to what we do for LSMR).
% 13 Jun 2013: (MZ) Added support for Method 22 for CME 338 project.
% 22 Nov 2013: (MAS) Various things polished up.  Method 22 made available.
% 23 Nov 2013: Restored options.backtrack.
%-----------------------------------------------------------------------

  PriLev   = options.Print;
  wait     = options.wait;
  if wait
    PriLev = 1;
  end

  if PriLev > 0
    fprintf('\n   --------------------------------------------------------')
    fprintf('\n   pdco.m                            Version of 23 Nov 2013')
    fprintf('\n   Primal-dual barrier method to minimize a convex function')
    fprintf('\n   subject to linear constraints Ax + r = b,  bl <= x <= bu')
    fprintf('\n                                                           ')
    fprintf('\n   Michael Saunders       SOL and ICME, Stanford University')
    fprintf('\n   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)')
    fprintf('\n                     Santiago Akle (ICME), Matt Zahr (ICME)')
    fprintf('\n   --------------------------------------------------------\n')
  end

  m = length(b);
  n = length(bl);


  %---------------------------------------------------------------------
  % Decode pdObj.
  %---------------------------------------------------------------------
  operator  =  isa(pdObj,'function_handle');
  explicitF = ~operator;       % Not a handle.  Must be a vector.
  
  if PriLev > 0
    if explicitF
      fprintf('\nThe objective is linear')
    else
      name = func2str(pdObj);
      fprintf('\nThe objective is defined by a function handle: ')
      if length(name) > 24, fprintf('\n   '), end
      fprintf('%s', name)
    end
  end

  %---------------------------------------------------------------------
  % Decode pdMat.
  %---------------------------------------------------------------------
  operator  =  isa(pdMat,'function_handle');
  explicitA = ~operator;       % Not a handle.  Must be a matrix.

  if PriLev > 0
    if explicitA
      nnzA   = nnz(pdMat);
      if issparse(pdMat)
        fprintf('\nThe matrix A is an explicit sparse matrix')
      else
        fprintf('\nThe matrix A is an explicit dense matrix' )
      end
      fprintf('\n\nm        = %8g     n        = %8g      nnz(A)  =%9g', m,n,nnzA)
    else
      name = func2str(pdMat);
      fprintf('\nThe matrix A is defined by a function handle: ')
      if length(name) > 24, fprintf('\n   '), end
      fprintf('%s', name)
      fprintf('\nm        = %8g     n        = %8g', m,n)
    end
  end

  normb  = norm(b ,inf);   normx0 = norm(x0,inf);
  normy0 = norm(y0,inf);   normz0 = norm(z0,inf);

  if PriLev > 0
    fprintf('\nmax |b | = %8g     max |x0| = %8.1e', normb , normx0)
    fprintf(                '      xsize   = %8.1e', xsize)
    fprintf('\nmax |y0| = %8g     max |z0| = %8.1e', normy0, normz0)
    fprintf(                '      zsize   = %8.1e', zsize)
  end

  %---------------------------------------------------------------------
  % Initialize.
  %---------------------------------------------------------------------
  true   = 1;
  false  = 0;
  zn     = zeros(n,1);
  nb     = n + m;
  CGitns = 0;
  inform = 0;
  % 07 Aug 2003: No need for next lines.
  %if length(d1)==1, d1 = d1*ones(n,1); end   % Allow scalar d1, d2
  %if length(d2)==1, d2 = d2*ones(m,1); end   % to mean d1*e, d2*e

  %---------------------------------------------------------------------
  % Grab input options.
  %---------------------------------------------------------------------
  maxitn    = options.MaxIter;
  featol    = options.FeaTol;
  opttol    = options.OptTol;
  steptol   = options.StepTol;
  stepSame  = options.StepSame;   % 1 means stepx==stepz
  x0min     = options.x0min;
  z0min     = options.z0min;
  mu0       = options.mu0;
  backtrack = options.backtrack;
  Method    = options.Method;     % 1=Cholesky  2=QR  3=LSMR  4=MINRES
                                  % 21=SQD(LU), 22=SQD(MA57)
  itnlim    = options.LSMRMaxIter * min(m,n);
  atol1     = options.LSMRatol1;  % Initial  atol
  atol2     = options.LSMRatol2;  % Smallest atol, unless atol1 is smaller
  conlim    = options.LSMRconlim;

  if Method==0
     if explicitA, Method = 1; else Method = 3; end
  end

  %---------------------------------------------------------------------
  % Set other parameters.
  %---------------------------------------------------------------------
  kminor    = 0;      % 1 stops after each iteration
  eta       = 1e-4;   % Linesearch tolerance for "sufficient descent"
  maxf      = 10;     % Linesearch backtrack limit (function evaluations)
  maxfail   = 1;      % Linesearch failure limit (consecutive iterations)
  bigcenter = 1e+3;   % mu is reduced if center < bigcenter

  % Parameters for LSMR and MINRES.
  atolmin   = eps;    % Smallest atol if linesearch back-tracks
  btol      = 0;      % Should be small (zero is ok)
  show      = false;  % Controls LSMR and MINRES iteration logs
  gamma     = max(d1);
  delta     = max(d2);

  if PriLev > 0
    fprintf('\n\nx0min    = %8g     featol   = %8.1e', x0min, featol)
    fprintf(                  '      d1max   = %8.1e', gamma)
    fprintf(  '\nz0min    = %8g     opttol   = %8.1e', z0min, opttol)
    fprintf(                  '      d2max   = %8.1e', delta)
    fprintf(  '\nmu0      = %8.1e     steptol  = %8g', mu0  , steptol)
    fprintf(                  '     bigcenter= %8g'  , bigcenter)

    fprintf('\n\nLSMR/MINRES:')
    fprintf('\natol1    = %8.1e     atol2    = %8.1e', atol1 , atol2 )
    fprintf(                  '      btol    = %8.1e', btol )
    fprintf('\nconlim   = %8.1e     itnlim   = %8g'  , conlim, itnlim)
    fprintf(                  '      show    = %8g'  , show )

    fprintf('\n\nMethod   = %8g     (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))\n', Method)

    if wait
      fprintf('\nReview parameters... then type "return"\n')
      keyboard
    end

    % if eta < 0
    %   fprintf('\n\nLinesearch disabled by eta < 0\n')
    % end
  end

  %---------------------------------------------------------------------
  % All parameters have now been set.
  % Check for valid Method.
  %---------------------------------------------------------------------
  time      = cputime;
  if PriLev > 0
    tic
  end
  
  if operator
    if Method==3 || Method ==4 %|| Method ==5
      % Relax
    else
      fprintf(['\n\n When A is an operator, we have to use Method = 3 or 4'])
      Method = 3;
    end
  end

  diagHess   = Method<=5;
  squareHess = Method==21 || Method==22;

  switch Method
  case 1
    solver  = '  Chol';  head3 = '     Chol';
  case 2
    solver  = '    QR';  head3 = '       QR';
  case 3
    solver  = '  LSMR';  head3 = '  atol   LSMR Inexact';
  case 4    
    solver  = 'MINRES';  head3 = '  atol MINRES Inexact'; 
% case 5    
%   solver  = '   PCG';  head3 = '  atol    PCG Inexact'; 
  case {21, 22}
    solver  = '   SQD';  head3 = '      SQD';
  otherwise
    error('Method must be 1, 2, 3, 4 or 21')
  end
    
  %---------------------------------------------------------------------
  % Categorize bounds and allow for fixed variables by modifying b.
  %---------------------------------------------------------------------
  [low,upp,fix,zlo,zup] = pdxxxbounds( bl,bu,PriLev );

  nfix = length(fix);
  if nfix > 0
    x1 = zn;   x1(fix) = bl(fix);
    r1 = pdxxxmat( pdMat, 1, m, n, x1 );
    b  = b - r1;
    % At some stage, might want to look at normfix = norm(r1,inf);
  end

  %---------------------------------------------------------------------
  % Scale the input data.
  % The scaled variables are
  %    xbar     = x/beta,
  %    ybar     = y/zeta,
  %    zbar     = z/zeta.
  % Define
  %    theta    = beta*zeta;
  % The scaled function is
  %    phibar   = ( 1   /theta) fbar(beta*xbar),
  %    gradient = (beta /theta) grad,
  %    Hessian  = (beta2/theta) hess.
  %---------------------------------------------------------------------
  beta   = xsize;   if beta==0, beta = 1; end    % beta scales b, x.
  zeta   = zsize;   if zeta==0, zeta = 1; end    % zeta scales y, z.
  theta  = beta*zeta;                            % theta scales obj.
  % (theta could be anything, but theta = beta*zeta makes
  % scaled grad = grad/zeta = 1 approximately if zeta is chosen right.)

  bl(fix)= bl(fix)/beta;
  bu(fix)= bu(fix)/beta;
  bl(low)= bl(low)/beta;
  bu(upp)= bu(upp)/beta;
  d1     = d1*( beta/sqrt(theta) );
  d2     = d2*( sqrt(theta)/beta );

  beta2  = beta^2;
  b      = b /beta;   y0     = y0/zeta;
  x0     = x0/beta;   z0     = z0/zeta;
 
  %---------------------------------------------------------------------
  % Initialize vectors that are not fully used if bounds are missing.
  %---------------------------------------------------------------------
  rL     = zn;   rU    = zn;
  cL     = zn;   cU    = zn;
  x1     = zn;   x2    = zn;
  z1     = zn;   z2    = zn;
  dx1    = zn;   dx2   = zn;
  dz1    = zn;   dz2   = zn;
  clear zn

  %---------------------------------------------------------------------
  % Initialize x, y, z1, z2, objective, etc.
  % 10 Aug 2003: z isn't needed here -- just at end for output.
  % 03 Jul 2008: Use pos to ensure that x = x1 for vanilla bounds.
  %              The linear constraints x(pos) - x1(pos) = bl(pos) = 0
  %              should then remain satisfied, so we'll have x(pos) > 0
  %              throughout.  At last, no danger of evaluating log(x)
  %              at negative x values.
  % 26 Feb 2010: zlo and zup now used in place of pos.
  %              See 26 Feb 2010 note above.
  %---------------------------------------------------------------------
  x      = x0;
  y      = y0;
  x(fix) = bl(fix);
  x(low) = max( x(low)          , bl(low));
  x(upp) = min( x(upp)          , bu(upp));
  x1(low)= max( x(low) - bl(low), x0min  );
  x2(upp)= max(bu(upp) -  x(upp), x0min  );
  z1(low)= max( z0(low)         , z0min  );
  z2(upp)= max(-z0(upp)         , z0min  );
  x(zlo) =   x1(zlo);
  x(zup) = - x2(zup);
  clear x0 y0 z0


  %%% "hess" can be a scalar, or a vector (diag(H)),
  %%% or an n x n sparse matrix.
  %%% If Method==21 or 22 (SQD), we need to make it the latter.
  
  if explicitF
     obj  = (pdObj'*x)*beta;  grad = pdObj;
     if diagHess
       hess = zeros(n,1);
     else
       hess = sparse(n,n);
     end
  else
    [obj,grad,hess] = pdObj(x*beta);     %#ok 
    [mH,nH] = size(hess);
    if diagHess
       if mH>1 && nH>1
          fprintf('\n Warning: Using only diagonal part of hess from pdObj\n')
          hess = diag(hess);
       end
    else
       if mH~=n && nH~=n
          error('For Method==21 || Method==22, pdObj must return a square (sparse) Hessian')
       end
    end
  end
  
  obj  = obj        /theta;               % Scaled obj.
  grad = grad*(beta /theta) + (d1.^2).*x; % grad includes x regularization.
  H    = hess*(beta2/theta);
  if diagHess
     H = H + (d1.^2);                     % H    includes x regularization.
  else
     H = H + sparse(1:n,1:n,(d1.^2),n,n); % H    includes x regularization.
  end

  %---------------------------------------------------------------------
  % Compute primal and dual residuals:
  %    r1 =  b - A*x - d2.^2*y
  %    r2 =  grad - A'*y + (z2-z1)
  %    rL =  bl - x + x1
  %    rU = -bu + x + x2
  %---------------------------------------------------------------------
  [r1,r2,rL,rU,Pinf,Dinf] = ...
      pdxxxresid1( pdMat,fix,low,upp, ...
                   b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 );

  %---------------------------------------------------------------------
  % Initialize mu and complementarity residuals:
  %    cL   = mu*e - X1*z1.
  %    cU   = mu*e - X2*z2.
  %
  % 25 Jan 2001: Now that b and obj are scaled (and hence x,y,z),
  %              we should be able to use mufirst = mu0 (absolute value).
  %              0.1 worked poorly on StarTest1 with x0min = z0min = 0.1.
  % 29 Jan 2001: We might as well use mu0 = x0min * z0min;
  %              so that most variables are centered after a warm start.
  % 29 Sep 2002: Use mufirst = mu0*(x0min * z0min),
  %              regarding mu0 as a scaling of the initial center.
  % 07 Aug 2003: mulast is controlled by opttol.
  %              mufirst should not be smaller.
  % 10 Aug 2003: Revert to mufirst = mu0 (absolute value).
  % 12 Jan 2006: If mu0 <= 0, reset it to the average complementarity.
  %---------------------------------------------------------------------
  if mu0 <= 0
    clow  = x1(low) .* z1(low);
    cupp  = x2(upp) .* z2(upp);
    mu0   = (sum(clow) + sum(cupp)) / (length(low) + length(upp));
    mu0   = mu0*0.1;
    clear clow cupp
  end

% mufirst = mu0*(x0min * z0min);
  mufirst = mu0;
  mulast  = 0.1 * opttol;
  mufirst = max( mufirst, mulast );
  mu      = mufirst;
  [cL,cU,center,Cinf,Cinf0] = ...
           pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
  fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );

  % Initialize other things.

  PDitns    = 0;
  converged = 0;
  atol      = atol1;
  atol2     = max( atol2, atolmin );
  atolmin   = atol2;
  btol      = atol;

  %  Iteration log.

  itncg   = 0;
  nfail   = 0;
  regterm = norm(d1.*x)^2  +  norm(d2.*y)^2;
  objreg  = obj  +  0.5*regterm;
  objtrue = objreg * theta;

  if PriLev > 0
    head1   = '\n\nItn   mu stepx stepz  Pinf  Dinf';
    head2   = '  Cinf   Objective    nf  center';
    fprintf([ head1 head2 head3 ])
    fprintf('\n%3g                 ', PDitns        )
    fprintf('%6.1f%6.1f' , log10(Pinf ), log10(Dinf))
    fprintf('%6.1f%15.7e', log10(Cinf0), objtrue    )
    fprintf('   %8.1f'   , center                   )

    if kminor
      fprintf('\n\nStart of first minor itn...\n')
      keyboard
    end
  end

  %---------------------------------------------------------------------
  % Main loop.
  %---------------------------------------------------------------------
  while ~converged
    PDitns = PDitns + 1;

    % 31 Jan 2001: Set atol according to progress, a la Inexact Newton.
    % 07 Feb 2001: 0.1 not small enough for Satellite problem.  Try 0.01.
    % 25 Apr 2001: 0.01 seems wasteful for Star problem.
    %              Now that starting conditions are better, go back to 0.1.

    r3norm = max([Pinf  Dinf  Cinf]);
    atol   = min([atol  r3norm*0.1]);
    atol   = max([atol  atolmin   ]);
    btol   = atol;

    if Method<=5  %%% diagHess
      %-----------------------------------------------------------------
      %  Solve (*) for dy.
      %-----------------------------------------------------------------
      %  Define a damped Newton iteration for solving f = 0,
      %  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
      %  to obtain the system
      %
      %  [-H2  A' ] [dx] = [w ],   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
      %  [ A  D2^2] [dy] = [r1]    w  = r2 - X1inv(cL + Z1 rL)
      %                                    + X2inv(cU + Z2 rU),
      %
      %  which is equivalent to the least-squares problem
      %
      %     min || [ D A']dy  -  [  D w   ] ||,   D = H2^{-1/2}.     (*)
      %         || [  D2 ]       [D2inv r1] ||
      %-----------------------------------------------------------------
      
      % For these methods to work, H must be diagonal 
      H(low)  = H(low) + z1(low)./x1(low);
      H(upp)  = H(upp) + z2(upp)./x2(upp);
      w       = r2;
      w(low)  = w(low) - (cL(low) + z1(low).*rL(low))./x1(low);
      w(upp)  = w(upp) + (cU(upp) + z2(upp).*rU(upp))./x2(upp);

      H      = 1./H;    % H is now Hinv (NOTE!)
      H(fix) = 0;
      D      = sqrt(H);
      pdDDD1 = D;
      

      if Method==1
        % --------------------------------------------------------------
        % Use chol to get dy
        % -------------------
        AD   = pdMat*sparse( 1:n, 1:n, D, n, n );
        ADDA = AD*AD' + sparse( 1:m, 1:m, (d2.^2), m, m );
        if PDitns==1, P = symamd(ADDA); end % Do ordering only once.

        [R,indef] = chol(ADDA(P,P));
        if indef
          fprintf('\n\n   chol says AD^2A'' is not pos def')
          fprintf('\n   Use bigger d2, or set options.Method = 2 or 3')
          inform = 4;
          break
        end

        % dy = ADDA \ rhs;
        rhs   = pdMat*(H.*w) + r1;
        dy    = R \ (R'\rhs(P));
        dy(P) = dy;
      
      elseif Method==2
        % --------------------------------------------------------------
        % Use QR to get dy
        % --------------------------------------------------------------
        DAt  = sparse( 1:n, 1:n, D, n, n )*(pdMat');
        if PDitns==1, P = colamd(DAt); end % Do ordering only once.

        if length(d2)==1
          D2 = d2*speye(m);
        else
          D2 = spdiags(d2,0,m,m);
        end

        DAt  = [ DAt;  D2     ];
        rhs  = [ D.*w; r1./d2 ];

        % dy = DAt \ rhs;
        [rhs,R] = qr(DAt(:,P),rhs,0);
        dy      = R \ rhs;
        dy(P)   = dy;

      elseif Method==3
        % --------------------------------------------------------------
        % Method=3.  Use LSMR (iterative solve) to get dy
        % --------------------------------------------------------------
        rhs    = [ D.*w; r1./d2 ];
        damp   = 0;

        if explicitA  % A is a sparse matrix.
                      % Construct diagonal preconditioner for LSMR
          precon = true;
          AD     = pdMat*diag(sparse(D));
          AD     = AD.^2;
          wD     = sum(AD,2); % Sum of sqrs of each row of AD.  %(Sparse)
          wD     = sqrt( full(wD) + (d2.^2) );                  %(Dense)
          pdDDD3 = 1./wD;
          clear AD wD

        else          % A is an operator. Can't use a preconditioner
          precon = false;  
          pdDDD3 = [];
        end

        mat_lsmr_handle = @(x,mode) ...
            pdxxxlsmrmat(mode, m, n, x, pdMat, Method, precon, pdDDD1, d2, pdDDD3);

        [dy, istop_lsmr, itncg, normr, normAr, normA, condA, normx] ...
            = lsmr(mat_lsmr_handle, rhs, damp, atol, btol, conlim,  ...
                   itnlim, [], show);  % dont set localSize

        %[ dy, istop, itncg, dr1norm, dr2norm, danorm, dacond, ...
        %  darnorm, dxnorm, dvar, dcov, outfo ] = ...
        %    lsqr( nb,m,matvechandle,rhs,damp,atol,btol,conlim, ...
        %          itnlim,show,[],[],stophandle);   % Don't want var or cov

        if precon, dy = pdDDD3 .* dy; end

        if istop_lsmr==3 || istop_lsmr==6 || istop_lsmr==7   % conlim or conda=1/eps or itnlim
          fprintf('\n    LSMR stopped early:  istop = %3d', istop_lsmr)
        end

        atolold   = atol;
        r3ratio   = normAr/fmerit;
        CGitns    = CGitns + itncg;
      
      elseif Method ==4
        % --------------------------------------------------------------
        % Method=4.  Use MINRES (iterative solve) to get dy.
        % --------------------------------------------------------------
        rhs    = pdxxxmat(pdMat,1,m,n,H.*w) + r1;
        damp   = 0;
        check  = 0;
        mat_minres_handle = @(x) pdxxxminresmat( m, n, x, pdMat, Method, H, d2 );

        if explicitA  % A is a sparse matrix.
                      % Construct diagonal preconditioner for MINRES
            precon = true;
            AD     = pdMat*diag(sparse(D));
            AD     = AD.^2;
            wD     = sum(AD,2);          % sparse
            wD     = full(wD) + d2.^2;   % dense
            pdDDD3 = 1./wD;  % Cast the vector into a diagonal operator
            pdDDD3 = @(x)(pdDDD3.*x);
        else
            precon = false;
            pdDDD3 = [];
        end

        [dy, istop_minres, itnm, normr, normAr, normA, condA, normx] = ...
            minres( mat_minres_handle, rhs, pdDDD3, damp, show, check, itnlim, atol );

        if istop_minres==-1 || istop_minres==5 || istop_minres==6   % conlim or itnlim
            fprintf('\n    MINRES stopped early:  istop = %3d', istop_minres)
        end

        atolold = atol;
        r3ratio = normr/fmerit; 
        CGitns  = CGitns + itnm;

      elseif Method==5
        % --------------------------------------------------------------
        % Method=5.  Use PCG (iterative solve) to get dy.
        % --------------------------------------------------------------
        error('PCG not implemented yet')

      end % computation of dy

      % dy is now known.  Get dx, dx1, dx2, dz1, dz2.

      grad      =   pdxxxmat( pdMat,2,m,n,dy );  % grad = A'dy
      grad(fix) =   0;    % Is this needed?      % grad is a work vector
      dx        =   H.*(grad - w);
      dx1(low)  = - rL(low) + dx(low);
      dx2(upp)  = - rU(upp) - dx(upp);
      dz1(low)  =  (cL(low) - z1(low).*dx1(low)) ./ x1(low);
      dz2(upp)  =  (cU(upp) - z2(upp).*dx2(upp)) ./ x2(upp);
    
      % See if atol and btol need to be reduced.

      if Method==3 || Method==4 || Method==5
         if atol > atolmin
            if r3ratio >= 0.001     % Accept dy but make next one more accurate.
               atol = atol*0.1;
            end
         end
      end
    end % if Method<=5


    if Method==21 || Method==22
      % ---------------------------------------------------------
      % Use SQD method to get dy and dx 
      % ---------------------------------------------------------
      H      = H + sparse(low,low, z1(low)./x1(low), n, n); 
      H      = H + sparse(upp,upp, z2(upp)./x2(upp), n, n); 
      w      = r2; 
      w(low) = w(low) - (cL(low) + z1(low).*rL(low))./x1(low);
      w(upp) = w(upp) + (cU(upp) + z2(upp).*rU(upp))./x2(upp);

      % Get rid of rows and columns of H corresponding to fixed
      % variables 
      
      if nfix > 0 
          [ih,jh,vh] = find(H); 
          clear H 
          for k=fix'
              vh(ih==k & ih~=jh) = 0;
              vh(jh==k & ih~=jh) = 0;
          end 
          H = sparse(ih,jh,vh); 
      end 
           
      if nfix==0
          K = [ -H     pdMat'
                pdMat  sparse(1:m, 1:m, d2.^2, m,m)];
      else
          if PDitns==1 
            Afree = pdMat;
            Afree(:,fix) = 0;
          end 
          K = [ -H     Afree'
                Afree  sparse(1:m, 1:m, d2.^2, m,m)];
      end

      if PDitns==1 && Method==21
          pSQD = symamd(K); 
      end % Do ordering only once.
        
      rhs = [w; r1]; 
      rhs(fix) = 0; 
     
      switch Method

        case 21
          % Use Matlab's old sparse LU (Gilbert and Peierls) on K. 
          % Note that K is symmetric quasi-definite (SQD).
          % If delta isn't too small, we can suppress row permutations
          % and still have a sufficiently stable method.

          thresh = eps;        % eps ~= 2e-16 suppresses partial pivoting
          [L,U,perm] = lu(K(pSQD,pSQD),thresh,'vector');
          if any(perm~=1:m+n)  % We trust this gives perm=1:n
              error('[L,U,perm]=lu(K(pSQD,pSQD),...) gave non-identity perm')
          end
          sqdsoln = U\(L\rhs(pSQD));
          sqdsoln(pSQD) = sqdsoln;

        case 22
          % Use MA57 via Matlab's sparse ldl.
          % MA57 uses multiple cores if available.
          thresh = 0;      % tells MA57 to keep its sparsity-preserving order
          [L,D,P,S] = ldl(K,thresh);
          if nnz(D)~=m+n
              error('[L,D,P,S] = ldl(K,0) gave non-diagonal D')
          end
          sqdsoln = S*(P*(L'\(D\(L\(P'*(S*rhs))))));
      end
      dx  = sqdsoln(1:n);
      dy  = sqdsoln(n+1:n+m); 
      dx1(low)  = - rL(low) + dx(low);
      dx2(upp)  = - rU(upp) - dx(upp);
      dz1(low)  =  (cL(low) - z1(low).*dx1(low)) ./ x1(low);
      dz2(upp)  =  (cU(upp) - z2(upp).*dx2(upp)) ./ x2(upp);
    end % if Method==21 || Method==22


    %-------------------------------------------------------------------
    % Find the maximum step.
    % 13 Aug 2003: We need stepxL, stepxU also to keep x feasible
    %              so that nonlinear functions are defined.
    % 18 Nov 2003: But this gives stepx = 0 for lptest.  (??)
    %--------------------------------------------------------------------
    stepx1 = pdxxxstep( x1(low), dx1(low) );
    stepx2 = pdxxxstep( x2(upp), dx2(upp) );
    stepz1 = pdxxxstep( z1(low), dz1(low) );
    stepz2 = pdxxxstep( z2(upp), dz2(upp) );
 %  stepxL = pdxxxstep(  x(low),  dx(low) );
 %  stepxU = pdxxxstep(  x(upp),  dx(upp) );
 %  stepx  = min( [stepx1,   stepx2,   stepxL,   stepxU] );
    stepx  = min( [stepx1,   stepx2] );
    stepz  = min( [stepz1,   stepz2] );
    stepx  = min( [steptol*stepx, 1] );
    stepz  = min( [steptol*stepz, 1] );
    if stepSame                      % For NLPs, force same step
      stepx = min( stepx, stepz );   % (true Newton method)
      stepz = stepx;
    end
   
    %-------------------------------------------------------------------
    % Backtracking linesearch.
    %-------------------------------------------------------------------
    fail     =  true;
    nf       =  0;

    while nf < maxf
      nf      = nf + 1;
      x1(low) = x1(low)  +  stepx * dx1(low);
      x2(upp) = x2(upp)  +  stepx * dx2(upp);
      z1(low) = z1(low)  +  stepz * dz1(low);
      z2(upp) = z2(upp)  +  stepz * dz2(upp);
      x       = x        +  stepx * dx;
      y       = y        +  stepz * dy;
      x(zlo)  =   x1(zlo);
      x(zup)  = - x2(zup);

      if explicitF
        obj = (pdObj'*x)*beta;  grad = pdObj;  % hess is already set
      else
        [obj,grad,hess] = pdObj(x*beta);       %#ok 
        if diagHess
           if mH>1 && nH>1
              hess = diag(hess);
           end
        end
      end
      obj        = obj /theta;
      grad       = grad*(beta /theta)  +  (d1.^2).*x;
      H          = hess*(beta2/theta);
      if diagHess
         H = H + (d1.^2);                     % H    includes x regularization.
      else
         H = H + sparse(1:n,1:n,(d1.^2),n,n); % H    includes x regularization.
      end

      [r1,r2,rL,rU,Pinf,Dinf] = ...
          pdxxxresid1( pdMat,fix,low,upp, ...
                       b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 );

      [cL,cU,center,Cinf,Cinf0] = ...
                  pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
      fmeritnew = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );
      step      = min( stepx, stepz );

      if ~backtrack | fmeritnew <= (1 - eta*step)*fmerit
         fail = false;
         break;
      end

      % Merit function didn't decrease.
      % Restore variables to previous values.
      % (This introduces a little error, but save lots of space.)
      
      x       = x        -  stepx * dx;
      y       = y        -  stepz * dy;
      x1(low) = x1(low)  -  stepx * dx1(low);
      x2(upp) = x2(upp)  -  stepx * dx2(upp);
      z1(low) = z1(low)  -  stepz * dz1(low);
      z2(upp) = z2(upp)  -  stepz * dz2(upp);
      x(zlo)  =   x1(zlo);
      x(zup)  = - x2(zup);

      % Back-track.
      % If it's the first time,
      % make stepx and stepz the same.

      if nf==1 && stepx~=stepz
         stepx = step;
      elseif nf < maxf
         stepx = stepx/2;
      end;
      stepz = stepx;
    end

    if fail
      fprintf('\n     Linesearch failed (nf too big)');
      nfail = nfail + 1;
    else
      nfail = 0;
    end

    %-------------------------------------------------------------------
    % Set convergence measures.
    %--------------------------------------------------------------------
    regterm = norm(d1.*x)^2  +  norm(d2.*y)^2;
    objreg  = obj  +  0.5*regterm;
    objtrue = objreg * theta;

    primalfeas    = Pinf  <=  featol;
    dualfeas      = Dinf  <=  featol;
    complementary = Cinf0 <=  opttol;
    enough        = PDitns>=       4;  % Prevent premature termination.
    converged     = primalfeas  &  dualfeas  &  complementary  &  enough;

    %-------------------------------------------------------------------
    % Iteration log.
    %-------------------------------------------------------------------
    if PriLev > 0
      str1    = sprintf('\n%3g%5.1f' , PDitns      , log10(mu)   );
      str2    = sprintf('%6.3f%6.3f' , stepx       , stepz       );
      if stepx < 0.001 || stepz < 0.001
        str2 = sprintf('%6.1e%6.1e' , stepx       , stepz       );
      end

      str3    = sprintf('%6.1f%6.1f' , log10(Pinf) , log10(Dinf));
      str4    = sprintf('%6.1f%15.7e', log10(Cinf0), objtrue     );
      str5    = sprintf('%3g%8.1f'   , nf          , center      );
      if center > 99999
         str5 = sprintf('%3g%8.1e'   , nf          , center      );
      end
      fprintf([str1 str2 str3 str4 str5])

      if     Method== 1
        if PDitns==1, fprintf(' %8g', nnz(R)); end
      elseif Method== 2
        if PDitns==1, fprintf(' %8g', nnz(R)); end
      elseif Method== 3
        fprintf(' %5.1f%7g%7.3f', log10(atolold), itncg, r3ratio)
      elseif Method== 4
        fprintf(' %5.1f%7g%7.3f', log10(atolold), itnm, r3ratio)
      elseif Method== 5
        fprintf(' %5.1f%7g%7.3f', log10(atolold), itnm, r3ratio)
 
      elseif Method==21
          if PDitns==1, fprintf(' %8g x2', nnz(L)); end
      elseif Method==22
          if PDitns==1, fprintf(' %8g', nnz(L)); end
      else
        %relax
      end

      %-----------------------------------------------------------------
      % Test for termination.
      %-----------------------------------------------------------------
      if kminor
        fprintf( '\nStart of next minor itn...\n')
        keyboard
      end
    end

    if converged
      if PriLev > 0
         fprintf('\n   Converged')
      end
    elseif PDitns >= maxitn
      if PriLev > 0
         fprintf('\n   Too many iterations')
      end
      inform = 1;
      break
    elseif nfail  >= maxfail
      if PriLev > 0
         fprintf('\n   Too many linesearch failures')
      end
      inform = 2;
      break
    elseif step   <= 1e-10
      fprintf('\nStep lengths too small')
      inform = 3;
      break
    else

      % Reduce mu, and reset certain residuals.

      stepmu  = min( stepx , stepz   );
      stepmu  = min( stepmu, steptol );
      muold   = mu;
      mumin   = 0.1*max([Pinf Dinf Cinf]);     % Target mu shouldn't be too small.
      mumin   = min(mu,mumin);                 % mu should be monotonic.
      
      mu      = mu - stepmu*mu;                % mu being treated as a variable.
      if center >= bigcenter, mu = muold; end  % Keep old mu if far from center.
      mu      = max(mu,mumin );                % 04 May 2008: No smaller than target.
      mu      = max(mu,mulast);                % 13 Jun 1998: No need for smaller mu.

      % mutrad = mu0*(sum(Xz)/n); % 24 May 1998: Traditional value, but
      % mu     = min(mu,mutrad ); % it seemed to decrease mu too much.

      mu      = max(mu,mulast);  % 13 Jun 1998: No need for smaller mu.
      [cL,cU,center,Cinf,Cinf0] = ...
               pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
      fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );

      % Reduce atol for LSMR (and MINRES).
      % NOW DONE AT TOP OF LOOP.

      atolold = atol;
      % if atol > atol2
      %   atolfac = (mu/mufirst)^0.25;
      %   atol    = max( atol*atolfac, atol2 );
      % end

      % atol = min( atol, mu );     % 22 Jan 2001: a la Inexact Newton.
      % atol = min( atol, 0.5*mu ); % 30 Jan 2001: A bit tighter

      % If the linesearch took more than one function (nf > 1),
      % we assume the search direction needed more accuracy
      % (though this may be true only for LPs).
      % 12 Jun 1998: Ask for more accuracy if nf > 2.
      % 24 Nov 2000: Also if the steps are small.
      % 30 Jan 2001: Small steps might be ok with warm start.
      % 06 Feb 2001: Not necessarily.  Reinstated tests in next line.


      if nf > 2  ||  step <= 0.01
        atol = atolold*0.1;
        btol = atol;
      end
    end
  end
  %---------------------------------------------------------------------
  % End of main loop.
  %---------------------------------------------------------------------

  % Reconstruct z.
  % Print statistics.
  
  x(fix) = 0;                 % Exclude x(fix) temporarily from |x|.
  z      = zeros(n,1);        % Exclude z(fix) also.
  z(low) = z1(low);
  z(upp) = z(upp) - z2(upp);

  if PriLev > 0
    fprintf('\n\nmax |x| =%10.3f', norm(x,inf))
    fprintf('    max |y| =%10.3f', norm(y,inf))
    fprintf('    max |z| =%10.3f', norm(z,inf))  % excludes z(fix)
    fprintf('   scaled')
  end

  bl(fix) = bl(fix)*beta;      % Unscale bl, bu, x, y, z.
  bu(fix) = bu(fix)*beta;
  bl(low) = bl(low)*beta;
  bu(upp) = bu(upp)*beta;

  x = x*beta;   y = y*zeta;   z = z*zeta;

  if PriLev > 0
    fprintf(  '\nmax |x| =%10.3f', norm(x,inf))
    fprintf('    max |y| =%10.3f', norm(y,inf))
    fprintf('    max |z| =%10.3f', norm(z,inf))  % excludes z(fix)
    fprintf(' unscaled')
  end

  % Reconstruct x(fix) and z(fix).
  r2 = pdxxxmat( pdMat,2,m,n,y );  % A'y
  if nfix > 0
    x(fix)  = bl(fix);
    z(fix)  = grad(fix) - r2(fix);        % z = grad - A'y
    z1(fix) = max( z(fix),0 );
    z2(fix) = max(-z(fix),0 );
  end

  % Reconstruct b.
  b = b*beta;
  if nfix > 0
    x1      = zeros(n,1);
    x1(fix) = bl(fix);
    r1      = pdxxxmat( pdMat,1,m,n,x1 );     % r1 = A*x1
    b       = b + r1;
    if PriLev > 0
       fprintf('\nmax |x| and max |z| exclude fixed variables')
    end
  end

  % Evaluate function at final point.
  % Recompute all of z.
  % 03 Apr 2010: z now includes (d1.^2).*x from the primal regularization.

  if explicitF
      grad = pdObj;     % This is all we need unless obj,grad,hess are output.
    % obj  = pdObj'*x;
    % grad = pdObj;
    % hess = zeros(n,1);
  else
      [obj,grad,hess] = pdObj( x );    %#ok
  end
  z    = grad + (d1.^2).*x - r2;     % z = grad + (d1.^2).*x - A'y

  time = cputime - time;
  if PriLev > 0
    fprintf('\nPDitns  =%10g %sitns =%10g    cputime =%10.1f', ...
               PDitns,solver,CGitns,time);
    pdxxxdistrib( abs(x),abs(z) );   % Private function
    toc
    if wait
      keyboard
    end
  end
%-----------------------------------------------------------------------
% End function pdco.m
%-----------------------------------------------------------------------


function [low,upp,fix,zlo,zup] = pdxxxbounds( bl,bu, PriLev )

% Categorize various types of bounds.
% pos overlaps with low.
% neg overlaps with upp.
% two overlaps with low and upp.
% fix and free are disjoint from all other sets.
%
% 26 Feb 2010: zlo and zup now point to variables with bl=0 or bu=0
%              respectively, with bl < bu.  Needed to keep x
%              strictly positive or negative.

  bigL = -9.9e+19;
  bigU =  9.9e+19;
  pos  =  find( bl==0    & bu>=bigU );
  neg  =  find( bl<=bigL & bu==0    );
  low  =  find( bl> bigL & bl< bu   );
  upp  =  find( bu< bigU & bl< bu   );
  two  =  find( bl> bigL & bu< bigU & bl< bu );
  fix  =  find( bl==bu );
  free =  find( bl<=bigL & bu>=bigU );
  zlo  =  find( bl==0    & bl< bu   );
  zup  =  find( bu==0    & bl< bu   );

  if PriLev > 0
    fprintf('\n\nBounds:\n  [0,inf]  [-inf,0]')
    fprintf('  Finite bl  Finite bu  Two bnds   Fixed    Free')
    fprintf('\n %8g %9g %10g %10g %9g %7g %7g', ...
        length(pos), length(neg), length(low),  ...
        length(upp), length(two), length(fix), length(free))
    fprintf('\n  [0, bu]  [bl,  0]  excluding fixed variables')
    fprintf('\n %8g %9g', length(zlo), length(zup))
  end
%-----------------------------------------------------------------------
% End private function pdxxxbounds
%-----------------------------------------------------------------------


function pdxxxdistrib( x,z )

% pdxxxdistrib(x) or pdxxxdistrib(x,z) prints the
% distribution of 1 or 2 vectors.
%
% 18 Dec 2000.  First version with 2 vectors.

  two  = nargin > 1;
  fprintf('\n\nDistribution of vector     x')
  if two, fprintf('         z'); end

  x1   = 10^(floor(log10(max(x)+eps)) + 1);
  z1   = 10^(floor(log10(max(z)+eps)) + 1);
  x1   = max(x1,z1);
  kmax = 10;

  for k = 1:kmax
    x2 = x1;    x1 = x1/10;
    if k==kmax, x1 = 0; end
    nx = length(find(x>=x1 & x<x2));
    fprintf('\n[%7.3g,%7.3g )%10g', x1, x2, nx);
    if two
      nz = length(find(z>=x1 & z<x2));
      fprintf('%10g', nz);
    end
  end

  disp(' ')
%-----------------------------------------------------------------------
% End private function pdxxxdistrib
%-----------------------------------------------------------------------


function y = pdxxxlsmrmat( mode, m, n, x, pdMat, Method, ...
                           precon, pdDDD1, d2, pdDDD3 )

% pdxxxlsmrmat is required by pdco.m when Method = 3.
% It forms Mx or M'x for some operator M that depends on Method below.
%
% m and n  are the dimensions of the LS problem that lsmr is solving
% (really n+m and m in terms of the m x n matrix A in pdco).
%
% pdMat is pdco's pdMat (defining A).
%
% rw contains parameters [explicitA Method LSdamp]
% from pdco.m to say which least-squares subproblem is being solved.
%
% pdDDD1 and pdDDD3 provide diagonal matrices for each value of Method.

%-----------------------------------------------------------------------
% 17 Mar 1998: First version to go with pdsco.m and lsqr.m.
% 01 Apr 1998: global pdDDD1 pdDDD3 now used for efficiency.
% 11 Feb 2000: Added diagonal preconditioning for LSQR, solving for dy.
% 14 Dec 2000: Added diagonal preconditioning for LSQR, solving for dx.
% 12 Feb 2001: Included in pdco.m as private function.
%              Specialized to allow only solving for dy.
% 03 Oct 2002: First version to go with pdco.m with general H2 and D2.
% 16 Oct 2002: pdMat is now the user's pdMat.
% 23 Apr 2008: (CMM) Added extra arguments Method, precon and
%              pdDDD1, pdDDD2, pdDDD3. pdDDD[123] are no longer
%              global variables. 
% 03 Apr 2012: Renamed to pdxxxlsmrmat.
%              Renamed nlsqr to n and mlsqr to m
%-----------------------------------------------------------------------

  if Method==3
    % The operator M is [ D1 A'; D2 ].
    if mode==1
      if precon, x = pdDDD3.*x; end
      t = pdxxxmat( pdMat, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
      y = [ (pdDDD1.*t); (d2.*x) ];

    else
      t = pdDDD1.*x(1:n);
      y = pdxxxmat( pdMat, 1, m, n, t );   % Ask 'aprod' to form y = A t.
      y = y + d2.*x(n+1:n+m);
      if precon, y = pdDDD3.*y; end
    end

  else
    error('Error in pdxxxlsmrmat: Only Method = 3 is allowed at present')
  end
%-----------------------------------------------------------------------
% End private function pdxxxlsmrmat
%-----------------------------------------------------------------------


function y = pdxxxminresmat( m, n, x, pdMat, Method, H, d2 )

% pdxxxminresmat is required by pdco.m when Method=4.
% It forms y=Mx or y=M'x for the operator M = A H A' + D2^2.
% m and n are the dimensions of the matrix A.
% pdMat is pdco's pdMat.

%-----------------------------------------------------------------------
% 04 Apr 2012: First version to go with pdco.m and minres.m.
%-----------------------------------------------------------------------

  if Method==4  % The operator M is [ A H A' + D2^2 ].
      t = pdxxxmat( pdMat, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
      t = t.*H;
      y = pdxxxmat( pdMat, 1, m, n, t );   %Ask 'aprod' to form y = At;
      y = y + (d2.^2).*x;
  else
      error('Error in pdxxxminresmat: Only Method 4 is allowed at present')
  end
%-----------------------------------------------------------------------
% End private function pdxxxminresmat
%-----------------------------------------------------------------------


function y = pdxxxmat( pdMat, mode, m, n, x )

%        y = pdxxxmat( pdMat, mode, m, n, x )
%    computes y = Ax (mode=1) or A'x (mode=2)
%    for a matrix A defined by pdco's input parameter pdMat.

%-----------------------------------------------------------------------
% 04 Apr 1998: Default A*x and A'*y function for pdco.m.
%              Assumed A was a global matrix pdAAA created by pdco.m
%              from the user's input parameter A.
% 16 Oct 2002: pdAAA eliminated to save storage.
%              User's parameter pdMat is now passed thru to here.
% 01 Nov 2002: Bug: feval had one too many parameters.
% 28 Nov 2003: y = pdMat'*x; replaced by y = (x'*pdMat)';
%-----------------------------------------------------------------------

  if isa(pdMat, 'function_handle')
    y = pdMat( mode, m, n, x );
  else
    if mode==1,  y = pdMat*x;  else  y = pdMat'*x;  end
  end
%-----------------------------------------------------------------------
% End private function pdxxxmat
%-----------------------------------------------------------------------


function fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU )

% Evaluate the merit function for Newton's method.
% It is the 2-norm of the three sets of residuals.

  f = [norm(r1)
       norm(r2)
       norm(rL(low))
       norm(rU(upp))
       norm(cL(low))
       norm(cU(upp))];
  fmerit = norm(f);
%-----------------------------------------------------------------------
% End private function pdxxxmerit
%-----------------------------------------------------------------------


function [r1,r2,rL,rU,Pinf,Dinf] =    ...
      pdxxxresid1( pdMat,fix,low,upp, ...
                   b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 )

% Form residuals for the primal and dual equations.
% rL, rU are output, but we input them as full vectors
% initialized (permanently) with any relevant zeros.
% 13 Aug 2003: z2-z1 coded more carefully
%              (although MATLAB was doing the right thing).
% 19 Nov 2003: r2(fix) = 0 has to be done after r2 = grad - r2;

  m       = length(b);
  n       = length(bl);
  x(fix)  = 0;
  r1      = pdxxxmat( pdMat, 1, m, n, x );
  r2      = pdxxxmat( pdMat, 2, m, n, y );

  r1      = b    - r1 - (d2.^2).*y;
  r2      = grad - r2;  % + (z2-z1);        % grad includes (d1.^2)*x
  r2(fix) = 0;
  r2(upp) = r2(upp) + z2(upp);
  r2(low) = r2(low) - z1(low);
  rL(low) = (  bl(low) - x(low)) + x1(low);
  rU(upp) = (- bu(upp) + x(upp)) + x2(upp);

  Pinf    = max([norm(r1,inf) norm(rL(low),inf) norm(rU(upp),inf)]);
  Dinf    =      norm(r2,inf);
  Pinf    = max( Pinf, 1e-99 );
  Dinf    = max( Dinf, 1e-99 );
%-----------------------------------------------------------------------
% End private function pdxxxresid1
%-----------------------------------------------------------------------


function [cL,cU,center,Cinf,Cinf0] = ...
      pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 )

% Form residuals for the complementarity equations.
% cL, cU are output, but we input them as full vectors
% initialized (permanently) with any relevant zeros.
% Cinf  is the complementarity residual for X1 z1 = mu e, etc.
% Cinf0 is the same for mu=0 (i.e., for the original problem).
% 12 Feb 2009: Beware: If all variables are free, Cinf0 = empty.
%                      Now changed to Cinf0 = 0;
% 03 Apr 2010: If all variables are free, we need to define center.
%              Set center = 1.0 arbitrarily.

  x1z1    = x1(low).*z1(low);
  x2z2    = x2(upp).*z2(upp);
  cL(low) = mu - x1z1;
  cU(upp) = mu - x2z2;

  maxXz   = max( [max(x1z1) max(x2z2)] );
  minXz   = min( [min(x1z1) min(x2z2)] );
  maxXz   = max( maxXz, 1e-99 );
  minXz   = max( minXz, 1e-99 );
  center  = maxXz / minXz;
  Cinf    = max([norm(cL(low),inf) norm(cU(upp),inf)]);
  Cinf0   = maxXz;
  if isempty(Cinf0)
    Cinf0  = 0;
    center = 1;
  end
%-----------------------------------------------------------------------
% End private function pdxxxresid2
%-----------------------------------------------------------------------

  
function step = pdxxxstep( x,dx )

% Assumes x > 0.
% Finds the maximum step such that x + step*dx >= 0.

  step     = 1e+20;
  blocking = find( dx < 0 );
  if ~isempty( blocking )
    steps  = x(blocking) ./ (- dx(blocking));
    step   = min( steps );
  end
%-----------------------------------------------------------------------
% End private function pdxxxstep
%-----------------------------------------------------------------------


%function [istop,atol,outfo] = pdxxxstop( istop,atol,arnorm,itn, ...
%                                         iatolmin,ir3norm )
%%-------------------------------------------------------------------
%% SPECIAL TEST THAT DEPENDS ON pdco.m.
%% pdMat in pdco   is  iw in lsqr.
%% dy              is  x
%% We allow for diagonal preconditioning in pdDDD3.
%%-------------------------------------------------------------------
%
%% This code gets input from lsqr (atol,arnorm,istop,itn)
%% and based on this data plus (iatolmin,ir3norm) adjusts (istop,atol).
%
%% 12 Feb 2001: atol can now be reduced and iterations continued
%%              if necessary.  outfo is a new
%%              problem-dependent parameter for such purposes.
%%              In this version they are specialized for pdco.m.
%% 23 Apr 2008: New function pdxxxstop handles adjusting atol inside lsqr.
%    
%
%  if istop > 0
%    r3new     = arnorm;
%    r3ratio   = r3new / ir3norm;
%    atolold   = atol;
%    atolnew   = atol;
%         
%    if atol > iatolmin
%      if     r3ratio <= 0.1     % dy seems good
%       % Relax
%      elseif r3ratio <= 0.5     % Accept dy but make next one more accurate.
%       atolnew = atolnew * 0.1;
%      else                      % Recompute dy more accurately
%       fprintf('\n                                ')
%       fprintf('                                ')
%       fprintf(' %5.1f%7g%7.3f', log10(atolold), itn, r3ratio)
%       atol    = atol * 0.1;
%       atolnew = atol;
%       istop   = 0;
%      end
%    end
%
%    outfo.atolold = atolold;
%    outfo.atolnew = atolnew;
%    outfo.r3ratio = r3ratio;
%  else
%    outfo = 0;
%  end
%%-----------------------------------------------------------------------
%% End private function pdxxxstop
%%-----------------------------------------------------------------------
