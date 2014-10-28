function pdcotestQP(m, n, alpha)
% m=50;  n=100;  pdcotestQP( m,n );
% Generates a random m by n QP problem
%   min x'Hx + c'x  st. Ax = b, bl < x < bu,
% and runs it on pdco.m.

%-----------------------------------------------------------------------
% 08 Oct 2002: Simple test program for pdco.m.
%              "A" is an explicit sparse matrix, not a function.
%              Michael Saunders, SOL, Stanford University.
% 16 Sep 2003: "c" is now an explicit vector.
%              It may be passed to pdco instead of function 'linobj'.
%-----------------------------------------------------------------------
  [A,b,bl,bu,c,d1,d2,H] = toydatatest(m,n,1); % Private function below
  %Add MOSEK path
  addpath /home/debasish/mosek/7/toolbox/r2013a/
  tic;
  mosekx = quadprog(H,c,[],[],A,b,bl,bu);
  mosekTime = toc;
  	
  %Add ECOS path
  addpath /home/debasish/ecos/matlab
  [ecosx, ~, ~, ~, ~, ecosTime] = ecosqp(H,c,[],[],A,b,bl,bu,ecosoptimset('verbose', 0, 'feastol', 1e-8));

  %Prepare gram matrix for equality constraint
  
  %Add ADMM path
  addpath /home/debasish/ecos/matlab/admm
  lambdaH = [H, A'; A, zeros(m)];

  rho = computeRho(lambdaH);
  fprintf('rho %g alpha %g\n', rho, alpha);
  
  tic;
  [xadmm history admmIters] = qpproximal(H, c, randn(1), A, b, bl, bu, rho, alpha, 0.0);
  admmTime = toc;
  
  tic;
  [xacc history accIters] = qpaccelerated(H, c, randn(1), A, b, bl, bu, rho, 0.0);
  accTime = toc;
  
% D  = sum(A,1);   D(find(D==0)) = 1;
% D  = sparse( 1:n, 1:n, 1./D, n, n );
% A  = A*D;                             % Normalize cols of A

%%%%%%% TEST OF FIXED VARIABLES ON pdcotestQP(5,10)
%if m==5 && n==10
%  bl(2) = 1.056179230657742;
%  bu(2) = bl(2);
%  bl(5) = 1.204010603155194;
%  bu(5) = bl(5);
%end

  options = pdcoSet;

  x0    = ones(n,1)/n;       % Initial x
  y0    = zeros(m,1);        % Initial y
  z0    = ones(n,1)/n;       % Initial z
  xsize = 5;                 % Estimate of norm(x,inf) at solution
  zsize = 8;                 % Estimate of norm(z,inf) at solution

  options.mu0       = 0e-0;  % An absolute value
  options.Method    = 21;    % SQD
  options.LSQRatol1 = 1e-8;  % For LPs, LSQR must solve quite accurately
  options.wait      = 0;     % Allow options to be reviewed before solve

  objectivefunction = @(x) deal(0.5*(x'*H*x) + c'*x, H*x + c, H);

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco(objectivefunction,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize );

  fprintf('mosek-ecos norm %g\n', norm(mosekx - ecosx, Inf));
  fprintf('mosek-pdco norm %g\n', norm(mosekx - x, Inf));
  fprintf('mosek-admm norm %g iters %g\n', norm(mosekx - xadmm, Inf), admmIters);
  fprintf('mosek-accadmm norm %g iters %g\n', norm(mosekx - xacc, Inf), accIters);
  fprintf('mosek %g ecos %g pdco %g admm %g admmacc %g\n', mosekTime, ecosTime, time, admmTime, accTime);
  
  %disp('Waiting in pdcotestQP')
  %keyboard                   % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% End function pdcotestQP
%-----------------------------------------------------------------------

function [A,b,bl,bu,c,d1,d2,H] = toydatatest(m,n, dense)

%        [A,b,bl,bu,c,d1,d2] = toydata( m,n );
%        defines an m by n matrix A, rhs vector b, and cost vector c
%        for use with pdco.m.

%-----------------------------------------------------------------------
% 12 Feb 2001: First version of toydata.m.
% 30 Sep 2002: pdco version modified for pdco.
% 16 Sep 2003: c now generated.  It may be passed to pdco.m.
% 23 Sep 2003: Simplified a bit to match pdcotestLS.m.
%-----------------------------------------------------------------------
  if dense > 0
     fprintf('Generating Dense Hessian and Inequalities\n');
  end
  
  %rand('state',10);

  density = 0.1;
  rc      = 1e-1;

  em      = ones(m,1);
  en      = ones(n,1);
  zn      = zeros(n,1);
  bigbnd  = 1e+30;

  if dense > 0
      A = rand(m,n);
  else                                 
      A = sprand(m,n,density,rc);
  end
  
  x       = en;

  gamma   = 1e-3;       % Primal regularization
  delta   = 1e-3;       % 1e-3 or 1e-4 for LP;  1 for Least squares.

  d1      = gamma;      % Can be scalar if D1 = d1*I.
  d2      = delta*em;

  b       = full(A*x);
  c       = rand(n,1);

  bl      = zn;         % x > 0
  bu      = 10*en;      % x < 10

  density = 0.05;

  if dense > 0
      H = rand(n, n);                                
  else                              
      H = sprand(n,n,density,rc);
  end
  
  H       = H'*H;
  H       = (H+H')/2;
% bl(1)   = - bigbnd;   % Test "free variable" (no bounds)
% bu(1)   =   bigbnd;

%-----------------------------------------------------------------------
% End private function toydata
%-----------------------------------------------------------------------
