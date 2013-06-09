function err = dimacs_errors(A,b,c,K,x,y,z)
% Computes DIMACS error measure. See the paper H.D. Mittelmann - "An
% independent benchmarking of SDP and SOCP solvers", Math. Program., 
% pp. 407-430, 2003.

if( nargin == 6 )
    z = c-A'*y;
end

lmin_x = min(x(K.f+(1:K.l)));
lmin_z = min(z(K.f+(1:K.l)));
for i = 1:length(K.q)
    coneidx = K.f + K.l + sum(K.q(1:(i-1))) + (1:K.q(i));
    lmin_x = min([lmin_x, x(coneidx(1))-norm(x(coneidx(2:end)),2)]);
    lmin_z = min([lmin_z, z(coneidx(1))-norm(z(coneidx(2:end)),2)]);
end

err(1) = norm(A*x-b) / (1 + norm(b,1));
err(2) = max([0, -lmin_x])/(1+norm(b,1));
err(3) = norm(A'*y + z - c)/norm(1+norm(c,1));
err(4) = max([0, -lmin_z])/(1+norm(c,1));
err(5) = (c'*x - b'*y)/(1 + abs(c'*x) + abs(b'*y));
err(6) = x'*z / (1 + abs(c'*x) + abs(b'*y));
