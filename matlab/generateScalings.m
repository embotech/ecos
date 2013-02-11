n = 4;
s0=1;
z0=1;
z1=1/n*randn(n-1,1); 
s1=1/n*randn(n-1,1);
s0 = norm(z1)+1e-6;
 
 s=[s0;s1]; 
 z=[z0;z1];
 [W,a,b,q,eta] = conelp_socscaling(s,z);
 
 w = q'*q;
 atilde = a^2 + w;
 c = 1 + a + w / (1+a);
 d = 1 + 2/(1+a) + w/(1+a)^2;
 
 V = eta^2 * [atilde, c*q'; c*q, eye(n-1) + d*(q*q')];
 Ve = eta^2 * [atilde, -c*q', 0; 
               c*q, eye(n-1), -sqrt(d)*q; 
               0, sqrt(d)*q', 1];
 Vprime = eta^2 * [atilde, c*q'; c*q, c*eye(n-1)];
 
 D = diag([a^2 + w; ones(n-1,1)]);
%  alpha = d;
%  beta = c^2/d;
%  v = [1; zeros(n-1,1)];
%  u = [c/d; q];
alpha = c^2/a^2/d;
u0 = a;
u1 = d/c*a;
beta = alpha*u0^2;
u = [u0; u1*q];
v = [1; zeros(n-1,1)];
Vnew = eta^2*(D + alpha*(u*u') - beta*(v*v'));


% search for other representation
u = sdpvar(n,1);
alpha = 1;%sdpvar(1);
% d = sdpvar(n,1);
%D = diag(d);
con =  [D alpha*u; alpha*u', -1] >= 0;
solvesdp(con,0)

 