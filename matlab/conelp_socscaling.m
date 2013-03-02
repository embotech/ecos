function [W,V,a,b,c,d,q,w,eta] = conelp_socscaling(s,z)
% Returns Nesterov-Todd scaling matrix for s,z in cone.

%% check whether s and z are in the cone.
s0 = s(1); s1 = s(2:end);
z0 = z(1); z1 = z(2:end);

sres = s0^2 - s1'*s1;
zres = z0^2 - z1'*z1;

assert( sres > 0, 's not in second-order cone');
assert( zres > 0, 'z not in second-order cone');


%% scalings
sbar = s./sqrt(sres);
zbar = z./sqrt(zres);
eta = (sres / zres)^(1/4);
gamma = sqrt( (1+sbar'*zbar)/2 );
wbar = 1/2/gamma*(sbar + [zbar(1); -zbar(2:end)]);
q = wbar(2:end);
a = wbar(1);
b = 1/(1+a);
W = eta*[a q'; q, eye(length(q)) + b*(q*q')];

w = q'*q;
c = 1 + a + w / (1+a);
d = 1 + 2/(1+a) + w/(1+a)^2;

V = eta^2*[a^2+w, c*q'; c*q, eye(length(sbar)-1)+d*(q*q')];