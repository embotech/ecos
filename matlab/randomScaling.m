function [Ve,W] = randomScaling(n)
% Generates a random scaling matrix (actually the square of it)
% with sparse expansion.
%
%   [Ve,V,W] = generateScalings(n) generates a random scaling matrix W, and
%   returns its square V = W^2, and the sparse expansion of V, Ve.
%


ok = 0;
while( ok == 0 )
try

z1=1/n*randn(n-1,1); 
z0=norm(z1) + 1;

s1=1/n*randn(n-1,1);
s0 = norm(s1)+1e-3;
 
s=[s0;s1]; 
z=[z0;z1];
%  [W,a,b,q,eta] = conelp_socscaling(s,z);

[W,V,a,b,c,d,q,w,eta] = conelp_socscaling(s,z);
catch
end
ok = 1;
end

%  w = q'*q;
%  atilde = a^2 + w;
%  c = 1 + a + w / (1+a);
%  d = 1 + 2/(1+a) + w/(1+a)^2;
%  
%  V = eta^2 * [atilde, c*q'; c*q, eye(n-1) + d*(q*q')];
%  Ve = eta^2 * [atilde, -c*q', 0; 
%                c*q, eye(n-1), -sqrt(d)*q; 
%                0, sqrt(d)*q', 1];
%  Vprime = eta^2 * [atilde, c*q'; c*q, c*eye(n-1)];
 
% D = diag([a^2 + w; ones(n-1,1)]);
%  alpha = d;
%  beta = c^2/d;
%  v = [1; zeros(n-1,1)];
%  u = [c/d; q];
%alpha = c^2/a^2/d;
%u0 = a;
%u1 = d/c*a;
%beta = alpha*u0^2;
%u = [u0; u1*q];
%v = [1; zeros(n-1,1)];

u0_lower = c^2/(1/w+d);
u0_upper = min([a^2+w, c^2/d]);
assert(u0_lower < u0_upper,'lower-upper thing does not hold');
u0 = sqrt( u0_lower + (u0_upper - u0_lower)/2);
u1 = c / u0;
v0 = 0;
v1 = sqrt(c^2/u0^2 - d);
assert(isreal(v1),'v1 is complex');
d1 = a^2 + w - u0^2;
assert(d1 > 0,'d1 is negative');

D = blkdiag(d1,eye(n-1));
u = [u0; u1*q];
v = [v0; v1*q];
Vnew = eta^2*(D + u*u' - v*v');
assert(all(eig(D-v*v'))>0,'D-vv'' has negative eigenvalues'); 

Ve = [D, v, u; v', 1, 0; u', 0, -1];
