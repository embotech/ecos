function [cscale,rscale] = gmscale(A,iprint,scltol)

%        [cscale,rscale] = gmscale(A,iprint,scltol);
%------------------------------------------------------------------
% gmscale (Geometric-Mean Scaling) finds the scale values for the
% m x n sparse matrix A.
%
% On entry:
% A(i,j)    contains entries of A.
% iprint    > 0 requests messages to the screen (0 means no output).
% scltol    should be in the range (0.0, 1.0).
%           Typically scltol = 0.9.  A bigger value like 0.99 asks
%           gmscale to work a little harder (more passes).
%
% On exit:
% cscale, rscale are column vectors of column and row scales such that
%           R(inverse) A C(inverse) should have entries near 1.0,
%           where R = diag(rscale), C = diag(cscale).
%
% Method:
% An iterative procedure based on geometric means is used,
% following a routine written by Robert Fourer, 1979.
% Several passes are made through the columns and rows of A.
% The main steps are:
%
%   1. Compute aratio = max_j (max_i Aij / min_i Aij).
%   2. Divide each row    i by  sqrt( max_j Aij * min_j Aij ).
%   3. Divide each column j by  sqrt( max_i Aij * min_i Aij ).
%   4. Compute sratio as in Step 1.
%   5. If  sratio < aratio * scltol,
%      set aratio = sratio and repeat from Step 2.
%
% To dampen the effect of very small elements, on each pass,
% a new row or column scale will not be smaller than sqrt(damp)
% times the largest (scaled) element in that row or column.
%
% Use of the scales:
% To apply the scales to a linear program,
%   min c'x  st  Ax = b,  l <= x <= u,
% we need to define "barred" quantities by the following relations:
%   A = R Abar C,  b = R bbar,  C cbar = c,
%   C l = lbar,    C u = ubar,  C x = xbar.
% This gives the scaled problem
%   min cbar'xbar  st  Abar xbar = bbar,  lbar <= xbar <= ubar.

% Maintainer:  Michael Saunders, Systems Optimization Laboratory,
%              Stanford University.
% 07 Jun 1996: First f77 version, based on MINOS 5.5 routine m2scal.
% 24 Apr 1998: Added final pass to make column norms = 1.
% 18 Nov 1999: Fixed up documentation.
% 26 Mar 2006: (Leo Tenenblat) First Matlab version based on Fortran version.
% 21 Mar 2008: (MAS) Inner loops j = 1:n optimized.
% 09 Apr 2008: (MAS) All loops replaced by sparse-matrix operations.
%              We can't find the biggest and smallest Aij
%              on each scaling pass, so no longer print them.
% 24 Apr 2008: (MAS, Kaustuv) Allow for empty rows and columns.
% 13 Nov 2009: gmscal.m renamed gmscale.m.
%------------------------------------------------------------------

  if iprint > 0
    fprintf('\ngmscale: Geometric-Mean scaling of matrix')
    fprintf('\n-------\n                 Max col ratio')
  end

  [m,n]   = size(A);
  A       = abs(A);    % Work with |Aij|
  maxpass = 10;
  aratio  = 1e+50;
  damp    = 1e-4;
  small   = 1e-8;
  rscale  = ones(m,1);
  cscale  = ones(n,1);

%---------------------------------------------------------------
% Main loop.
%---------------------------------------------------------------
  for npass = 0:maxpass

    % Find the largest column ratio.
    % Also set new column scales (except on pass 0).

    rscale(rscale==0) = 1;
    Rinv    = diag(sparse(1./rscale));
    SA      = Rinv*A;
    [I,J,V] = find(SA);
    invSA   = sparse(I,J,1./V,m,n);
    cmax    = full(max(   SA))';   % column vector
    cmin    = full(max(invSA))';   % column vector
    cmin    = 1./(cmin + eps);
    sratio  = max( cmax./cmin );   % Max col ratio
    if npass > 0
      cscale = sqrt( max(cmin, damp*cmax) .* cmax );
    end

    if iprint > 0
      fprintf('\n  After %2g %19.2f', npass, sratio)
    end

    if npass >= 2 && sratio >= aratio*scltol, break; end
    if npass == maxpass, break; end
    aratio  = sratio;

    % Set new row scales for the next pass.

    cscale(cscale==0) = 1;
    Cinv    = diag(sparse(1./cscale));
    SA      = A*Cinv;                  % Scaled A
    [I,J,V] = find(SA);
    invSA   = sparse(I,J,1./V,m,n);
    rmax    = full(max(   SA,[],2));   % column vector
    rmin    = full(max(invSA,[],2));   % column vector
    rmin    = 1./(rmin + eps);
    rscale  = sqrt( max(rmin, damp*rmax) .* rmax );
  end
%---------------------------------------------------------------
% End of main loop.
%---------------------------------------------------------------

% Reset column scales so the biggest element
% in each scaled column will be 1.
% Again, allow for empty rows and columns.

  rscale(rscale==0) = 1;
  Rinv    = diag(sparse(1./rscale));
  SA      = Rinv*A;
  [I,J,V] = find(SA);
  cscale  = full(max(SA))';   % column vector
  cscale(cscale==0) = 1;

% Find the min and max scales.

  if iprint>0
    [rmin,imin] = min(rscale);
    [rmax,imax] = max(rscale);
    [cmin,jmin] = min(cscale);
    [cmax,jmax] = max(cscale);

    fprintf('\n\n  Min scale               Max scale')
    fprintf('\n  Row %6g %9.1e    Row %6g %9.1e', imin, rmin, imax, rmax)
    fprintf('\n  Col %6g %9.1e    Col %6g %9.1e', jmin, cmin, jmax, cmax)
  end

% end of gmscale
