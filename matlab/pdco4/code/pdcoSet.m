function options = pdcoSet(opts)

% pdcoSet creates or alters an options structure for pdco.m.
%
%   options = pdcoSet; (with no input arguments)
%   creates a structure with all fields set to their default values.
%   Each field is an option (also called a parameter).
%
%   pdcoSet (with no input or output arguments)
%   displays all options and their default values.
%
%   options = pdcoSet(opts); 
%   creates a structure with all fields set to their default values,
%   except valid fields in the structure "opts" replace the defaults.
%
% options.MaxIter     Maximum iterations of the primal-dual barrier method.
%                     Most problems should solve within 30 PDitns.
% options.FeaTol      Accuracy for satisfying Ax + D2 r = b, A'y + z = gobj
%                     and x - x1 = bl, x + x2 = bu, where x1, x2 > 0.
%                     1e-6 is typically small enough.
%                     1e-5 may be acceptable also.
% options.OptTol      Accuracy for satisfying x1.*z1 = 0, x2.*z2 = 0,
%                     where z = z1 - z2 and z1, z2 > 0.
%                     Typically the same as Featol.
% options.Print       1 gives output (default).  0 suppresses output.
% options.StepTol     (between 0 and 1): Controls how close each an
%                     x or z may be to reaching a bound at each step.
%                     For safety, should not be bigger than 0.99 (say)
%                     for nonlinear problems.
% options.StepSame    1 (true) if stepx and stepz should be the same
%                     (gives true Newton method for nonlinear problems);
%                     0 (false) if stepx and stepz may be different
%                     (gives better performance for linear programs).
% options.x0min       Min distance between x0 and bl or bu  AFTER SCALING.
%                     1.0 is about right for cold starts.
%                     0.1 or 0.01 may be ok if x0, z0 are "good".
% options.z0min       Min distance between abs(z0) and zero AFTER SCALING,
%                     when z0 is used to initialize z1 > 0 and z2 > 0.
%                     Typically the same as x0min.
% options.mu0         Initial mu (ABSOLUTE VALUE) for solving scaled problem.
% options.backtrack   Specifies whether to do a backtracking linesearch
%                     on the norm of the residuals for the Newton system.
% options.Method      Specifies how each search direction (dx,dy,dz1,dz2)
%                     should be computed.  Several methods exist for
%                     experimental purposes.  If A has fewer rows than columns
%                     (m < n) we usually solve for dy first (most of the work)
%                     and then get the other components cheaply.
%
%  Method  Solve for  Using                                         Implemented?
%     0       dy      Method=1 or Method=3 if A = matrix/operator resp. Yes
%     1       dy      Sparse Cholesky on (A D^2 A' + D2^2 I).           Yes
%     2       dy      Sparse QR on corresponding least-squares problem  Yes
%     3       dy      LSMR on least-squares problem                     Yes
%     4       dy      MINRES                                            Yes
%
%    11       dx      Sparse Cholesky on (D A'A D  + D2^2 I).           No
%    12       dx      Sparse QR on corresponding least-squares problem  No
%    13       dx      LSQR on least-squares problem                     No
%
%    21    dx,dy      Sparse LU on 2x2 KKT-type system                  Yes
%    22    dx,dy      Sparse indefinite LDL' on 2x2 KKT-type system     Yes
%    23    dx,dy      SYMMLQ on same system                             No
%
%    31    dx,dy,dz1  Sparse LU on 3x3 system (only for problems with   No
%                     with vanilla bounds: 0 < x < inf).
%                     This is a HUGE system, but it is relatively
%                     well-conditioned compared to the other systems.
%
%    41 dx,dy,dz1,dz2 Sparse LU on 4x4 system with general bounds.    Deleted
%                     This is an even HUGER system.
%
%    If A is an explicit sparse matrix, all methods are applicable.
%     1 is usually best (e.g. for LPs).
%     2 may be more reliable; it's there for checking.
%     3 is sometimes more efficient (e.g. for entropy problems).
%       Diagonal preconditioning is possible.
%     4 may be more efficient under the same circumstances.
%
%    If A is an operator,
%     3 or 4 must be used.  Diagonal preconditioning is not possible.
%
%    Notes:
%    Method =  1      On ENTROPY.big, symamd can find an ordering
%                     but chol never finishes.
%
% The following options control LSMR when Method = 3:
%
% options.LSMRMaxIter * min(m,n) is the maximum LSMR  iterations.
% options.LSMRatol1   is the starting value of the LSMR accuracy 
%                     tolerance "atol" (if LSmethod = 3).
%                     1e-3 or 1e-4 sometimes works.
%                     1e-8 may be needed for LPs.
%                     In general, if max(Pinf,Dinf,Cinf) doesn't decrease
%                     every iteration, set atol1 to a smaller value.
% options.LSMRatol2   is the smallest value atol is reduced to.
% options.LSMRconlim  shuts LSMR down early if its matrix is ill-conditioned.
%
% options.wait = 0    means pdco should proceed to solve the problem.
%              = 1    means pdco should pause to allow interactive resetting
%                     of some of the parameters.


% 09 Dec 2008: pdcoSet follows Michael Friedlander's as_setparms.m.
% 02 Jun 2011: options.backtrack added.  The default value 0
%              turns the search OFF because it seems to inhibit
%              convergence on some entropy problems arising in a
%              fixed-point iteration in systems biology.
% 03 Apr 2012: Changed LSQR parameters to LSMR parameters; added MINRES.
% 05 Apr 2012: Added default Method = 0.

% Set default options.
  options.MaxIter      =    30;
  options.FeaTol       =  1e-6;
  options.OptTol       =  1e-6;
  options.Print        =     1;
  options.StepTol      =  0.99;
  options.StepSame     =     1;  % 1 for stepx == stepz (NLPs)
  options.x0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
  options.z0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
  options.mu0          =  1e-1;  % < 1.0 better than 1.0?
  options.backtrack    =     0;  % 0 disables the backtracking linesearch
  options.Method       =     0;  % Will change to 1 or 3
  options.LSMRMaxIter  =  10.0;
  options.LSMRatol1    = 1e-10;
  options.LSMRatol2    = 1e-15;  % 
  options.LSMRconlim   = 1e+12;  % Somewhere between e+8 and e+16
  options.wait         =     0;
  options.NOTE         = 'LSMRMaxIter is scaled by the matrix dimension';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
  if nargin == 0 || isempty(opts)
    if nargout == 0    % Display options.
      disp('pdco default options:')
      disp( options )
    end
    return
  end

% List of valid field names
  vfields = fieldnames( options );

% Grab valid fields from user's opts
  for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
      options.(field) = opts.(field);
    end
  end
