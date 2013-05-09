function [x,y,z,nitref] = conelp_solve(L,D,P,PL,QL, bx,by,bz, A,G,V, dims, nItref, LINSOLVER, LINSYSACC) %A,G,Gtilde,V,dims, K, Dreg, nItref, ppp)%,W,G,scaling,Gtilde,Winv,A,V)
% Solve KKT system using the cached factors L,D and permutation matrix P.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2013.

n = length(bx);
p = length(by);

% determine how many new rows are to be added
switch( LINSOLVER )
    case 'backslash',    Nstretch = 0;
    case 'rank1updates', Nstretch = 1;
    case 'ldlsparse',    Nstretch = 2;
end

% prepare RHS with zeros at appropriate places
bztilde = conelp_stretch(bz, dims, Nstretch);
RHS = [bx; by; bztilde];
PRHS = full(RHS); PRHS = PRHS(P);

% initial solve
switch( LINSOLVER )
    case 'backslash'
        dx(P,1) = L'\(D\(L\PRHS));
    
    case 'rank1updates'
        dx(P,1) = conelp_lowranksolve(L,D,PL,QL,PRHS);
        
    case 'ldlsparse'
        u = conelp_forwardsub(L,PRHS);
        v = conelp_byDiag(D,u);
        dx(P,1) = conelp_backwardsub(L',v);
        
    otherwise, error('Unknown linear solver');
end

% iterative refinement
bnorm = 1+norm(PRHS,inf);
for i = 1:nItref
    
    % variables
    x = dx(1:n);
    y = dx(n+1:n+p);
    ztilde = dx(n+p+1:end);
    z = conelp_unstretch(ztilde,dims, Nstretch);
    
    % errors
    ex = bx - A'*y - G'*z;
    ey = by - A*x;
    eztilde = conelp_stretch(bz - G*x,dims,2) + V*ztilde;
%     e = [ex; ey; conelp_stretch(ez,dims,Nstretch)];
    e = [ex; ey; eztilde];
%         fprintf('||ex||=%4.2e  ||ey||=%4.2e  ||ez||=%4.2e  (k=%d)\n', norm(ex)/bnorm, norm(ey)/bnorm, norm(eztilde)/bnorm, i);
        if(norm(ex,inf)/bnorm < LINSYSACC && ...
           norm(ey,inf)/bnorm < LINSYSACC && ...
           norm(eztilde,inf)/bnorm < LINSYSACC ), break; end

%         fprintf('*');
    
    % solve for correction
    switch( LINSOLVER )
        case 'backslash'
            dx(P,1) = dx(P,1) + L'\(D\(L\e(P)));
        
        case 'rank1updates'
            dx(P,1) = dx(P,1) + conelp_lowranksolve(L,D,PL,QL,e(P));
            
        case 'ldlsparse'
            u = conelp_forwardsub(L,e(P));
            v = conelp_byDiag(D,u);
            dx(P,1) = dx(P,1) + conelp_backwardsub(L',v);
            
        otherwise, error('Unknown linear solver');
    end
end
% fprintf('---\n');
% copy out variables
x = dx(1:n);
y = dx(n+1:n+p);
z = conelp_unstretch(dx(n+p+1:end),dims,Nstretch);
nitref = i;
% fprintf('\n');

