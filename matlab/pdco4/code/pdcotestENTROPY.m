function [x,y,z,inform] = pdcotestENTROPY(fname)
%        [x,y,z,inform] = pdcotestENTROPY(fname);
%        [x,y,z,inform] = pdcotestENTROPY('ENTROPY.small');
% Loads network data from file fname
% and solves a network problem with entropy objective,
%    min sum xj*ln(xj)  s.t.  Ax = b, x > 0.
%
% The network data is a file of (i,j) pairs.
% The corresponding variables xij must satisfy
%    sum(i) xik  -  sum(j) xkj  = 0,
%                  sum(ij) xij  = 1.

%-----------------------------------------------------------------------
% 08 Feb 2002: First version of ent3.m
%              derived from John Tomlin's AMPL model.
%              Network data built in as an explicit matrix.
%              Michael Saunders (SOL) and John Tomlin (IBM).
% 09 Feb 2002: Load (i,j) data from a file and use it to create A.
% 13 Mar 2002: Expand A to ensure feasibility.
% 28 May 2003: First version of ent5.m -- calls pdco instead of pdsco.
% 13 Aug 2003: mu0 is now an absolute number.
% 17 Sep 2003: pdcotestENTROPY.m made from ent5.m.
% 19 Apr 2010: The entropy objective is now passed as a function handle
%              to pdco version 03 Apr 2010.
%-----------------------------------------------------------------------

A = load(fname);      % A(:,1) -> +1,   A(:,2) -> -1.
k = min(A(:,1));
m = max(A(:,1));
if k==0   % Adjust for 0 base.
   A = A + 1;   m = m+1;
end
n = size(A,1);        % n = number of links
I = (1:m)';           % m = number of nodes
J = (1:n)';           % Column indices
m1= m+1;              % Include extra constraint, sum(x) = 1
m2= m+2;              % Row of -1s under  I
m3= m+3;              % Row of  1s under -I
em= ones(m,1);
en= ones(n,1);
nA= n + 2*m;   enA = ones(nA,1);

A = sparse([A(:,1); A(:,2);   I;     I;  m1*em;  m2*em;  m3*enA], ...
           [  J   ;   J   ; n+I; n+m+I;    n+I;  n+m+I; (1:nA)'], ...
           [  en  ;  -en  ;  em;   -em;    -em;     em;   enA  ]);

% For ent2.dat, A will now be this:
% [ 1  1  0  0  0  0  0  0  0  0  0  0  0  0 -1 1          -1
%  -1  0  1  1  0  0  0  0  0  0  0  0  0  0  0   1        -1
%     -1  0  0  1  1  0  0  0  0  0  0  0  0  0     .         .
%        -1  0 -1  0  1  1  0  0  0  0  0  0  0      .         .
%           -1  0 -1  0  0  1  1  1  0  0  0  0       .         .
%                    -1  0 -1  0  0  1  0  0  0
%                       -1  0 -1  0  0  1  0  0
%                                -1  0  0  1  0
%                                   -1 -1 -1  1           1        -1
%                                               -1 -1 ...-1
%                                                           1 1 ... 1
%   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 ..  1 1 1 ... 1]

alpha = 0.9;
b     = [zeros(m,1); -(1-alpha); (1-alpha); 1];
n     = n + 2*m;
m     = m3;
en    = ones(n,1);

bl    = zeros(n,1);
bu    = en*1e+20;

xsize = 5/n;               % A few elements of x are much bigger than 1/n.
xsize = min(xsize,1);      % Safeguard for tiny problems.
zsize = 1;                 % This makes y (sic) about the right size.
			   % 10 makes ||y|| even closer to 1,
			   % but for some reason doesn't help.

x0min = xsize;             % Applies to scaled x1, x2
z0min = zsize;             % Applies to scaled z1, z2

x0    = en*xsize;          % 
y0    = zeros(m,1);
z0    = en*z0min;          % z is nominally zero (but converges to mu/x)

d1    = 0;                 % gamma. 1e-3 is normal.  0 seems fine for entropy
d2    = 1e-3;              % delta

options = pdcoSet;
options.MaxIter      =    50;
options.FeaTol       =  1e-6;
options.OptTol       =  1e-6;
options.x0min        = x0min;  % This applies to scaled x1, x2.
options.z0min        = z0min;  % This applies to scaled z1, z2.
options.mu0          =  1e-5;  % 09 Dec 2005: BEWARE: mu0 = 1e-5 happens
                               %    to be ok for the entropy problem,
                               %    but mu0 = 1e-0 is SAFER IN GENERAL.

options.Method       =     3;  % 1=Chol  2=QR  3=LSQR
options.LSMRatol1    =  1e-3;
options.LSMRatol2    =  1e-6;
options.wait         =     1;

[x,y,z,inform] = pdco(@entropy,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize);

%-----------------------------------------------------------------------
% End function pdcotestENTROPY
%-----------------------------------------------------------------------
