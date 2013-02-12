function [x,y,z] = conelp_solve(L,D,P, bx,by,bz, A,G,V, dims, nItref) %, LINSOLVER, A,G,Gtilde,V,dims, K, Dreg, nItref, ppp)%,W,G,scaling,Gtilde,Winv,A,V)
% Solve KKT system using the cached factors L,D and permutation matrix P.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2013.

n = length(bx);
p = length(by);

% prepare RHS with zeros at appropriate places
bztilde = conelp_stretch(bz, dims);
RHS = [bx; by; bztilde];
PRHS = full(RHS); PRHS = PRHS(P);

% solve
u = conelp_forwardsub(L,PRHS);
v = conelp_byDiag(D,u);
dx(P,1) = conelp_backwardsub(L',v);

% iterative refinement
% bnorm = norm(PRHS);
for i = 1:nItref
    
    % variables
    x = dx(1:n);
    y = dx(n+1:n+p);
    z = conelp_unstretch(dx(n+p+1:end),dims);
%     z = dx(n+p+1:end);
    
    % errors
    ex = bx - A'*y - G'*z;
    ey = by - A*x;
    ez = bz - G*x + V*z;
    e = [ex; ey; conelp_stretch(ez,dims)];
%     fprintf('||ex||=%4.2e  ||ey||=%4.2e  ||ez||=%4.2e  (k=%d)\n', norm(ex), norm(ey), norm(ez), i);
%     if(norm(ex)/bnorm < 1e-9 && norm(ey)/bnorm < 1e-9 && norm(ez)/bnorm < 1e-13 ), break; end    
%     fprintf('*');
    % solve
    u = conelp_forwardsub(L,e(P));
    v = conelp_byDiag(D,u);
    dx(P,1) = dx(P,1) + conelp_backwardsub(L',v);
end

% copy out variables
x = dx(1:n);
y = dx(n+1:n+p);
% z = dx(n+p+1:end);
z = conelp_unstretch(dx(n+p+1:end),dims);

% fprintf('\n');

