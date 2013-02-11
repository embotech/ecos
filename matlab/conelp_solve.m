function [x,y,z] = conelp_solve(L,D,P, bx,by,bz, A,Gtilde,V, dims, nItref) %, LINSOLVER, A,G,Gtilde,V,dims, K, Dreg, nItref, ppp)%,W,G,scaling,Gtilde,Winv,A,V)
% Solve KKT system using the cached factors L,D and permutation matrix P.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2013.

n = length(bx);
p = length(by);

% prepare RHS with zeros at appropriate places
% bztilde = conelp_stretch(bz, dims);
RHS = [bx; by; bz];
PRHS = full(RHS); PRHS = PRHS(P);

% solve
u = conelp_forwardsub(L,PRHS);
v = conelp_byDiag(D,u);
dx(P,1) = conelp_backwardsub(L',v);

% iterative refinement
bnorm = norm(PRHS);
for i = 1:nItref
    
    % variables
    x = dx(1:n);
    y = dx(n+1:n+p);
    z = dx(n+p+1:end);
    
    % errors
    ex = bx - A'*y - Gtilde'*z;
    ey = by - A*x;
    ez = bz - Gtilde*x + V*z;
    e = [ex; ey; ez];
%     fprintf('||ex||=%4.2e  ||ey||=%4.2e  ||ez||=%4.2e  (k=%d)\n', norm(ex), norm(ey), norm(ez), i);
    if(norm(ex)/bnorm < 1e-11 && norm(ey)/bnorm < 1e-11 && norm(ez)/bnorm < 1e-11 ), break; end    
    fprintf('*');
    % solve
    u = conelp_forwardsub(L,e(P));
    v = conelp_byDiag(D,u);
    dx(P,1) = dx(P,1) + conelp_backwardsub(L',v);
end

% copy out variables
x = dx(1:n);
y = dx(n+1:n+p);
z = dx(n+p+1:end);
% z = conelp_unstretch(dx(n+p+1:end),dims);

fprintf('\n');

%% DEPRECATED
% if( nargin < 14 )
%     Dreg = eye(size(V));
% end
%
% if( nargin < 15 )
%     nItref = 3;
% end
%
% nItref = 3;
%
% switch( lower(LINSOLVER) )
%
%     case 'backslash'
%         % prepare RHS with zeros at appropriate places
%         RHS = [bx; by; bz(1:dims.l)];
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + (1:dims.q(k));
%             RHS = [RHS; bz(coneidx); 0];
%         end
%
%         % solve
%         u = L\(P'*RHS); v = D\u; w = L'\v; dx = P*w;
%
%         % copy out variables
%         x = dx(1:n); y = dx(n+1:n+p);
%         z = dx(n+p+1:n+p+dims.l);
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + k-1 + (1:dims.q(k));
%             z = [z; dx(n+p+coneidx)];
%         end
%
%     case 'backslash2'
%         RHS = [bx + G'/W/W*bz; by];
%         u = L\(P'*RHS); v = D\u; w = L'\v; dx = P*w;
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         z = W'\(W\G*x - W\bz);
%
%     case 'backslash3'
%         %         RHS = [bx; by; W\bz];
%         RHS = [bx; by; conelp_byW(scaling,bz,dims)];
%         u = L\(P'*RHS); v = D\u; w = L'\v; dx = P*w;
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         %         z = W\dx(n+p+1:end);
%         z =  conelp_byW(scaling,dx(n+p+1:end),dims);
%
%     case 'ldlsparse'
%
%         % prepare RHS with zeros at appropriate places
%         bztilde = bz(1:dims.l);
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + (1:dims.q(k));
%             bztilde = [bztilde; bz(coneidx); 0];
%         end
%         RHS = [bx; by; Dreg\bztilde];
%         PRHS = full(RHS); PRHS = PRHS(P);
%
%         % solve
%         u = conelp_forwardsub(L,PRHS);
%         v = conelp_byDiag(D,u);
%         dx(P,1) = conelp_backwardsub(L',v);
%
%         % iterative refinement
%         for i = 1:nItref
%
%                 % variables
%                 x = dx(1:n);
%                 y = dx(n+1:n+p);
%                 Dz = dx(n+p+1:end);
%                 z = Dreg\Dz;
%
%                 % errors
%                 ex = bx - A'*y - Gtilde'*z;
%                 ey = by - A*x;
%                 ez = Dreg\(bztilde - Gtilde*x + V*z);
%                 e = [ex; ey; ez];
%
% %                 fprintf('||ex||=%4.2e  ||ey||=%4.2e  ||ez||=%4.2e  (k=%d)\n', norm(ex), norm(ey), norm(ez), i);
%
%                 % solve
%                 u = conelp_forwardsub(L,e(P));
%                 v = conelp_byDiag(D,u);
%                 dx(P,1) = dx(P,1) + conelp_backwardsub(L',v);
%         end
%
%
%         % copy out variables
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         dx(n+p+1:end) = Dreg\dx(n+p+1:end); %  z = D \ D*z
%         z = dx(n+p+1:n+p+dims.l);
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + k-1 + (1:dims.q(k));
%             z = [z; dx(n+p+coneidx)];
%         end
%
%
%     case 'ldlsparse4'
%         % prepare RHS with zeros at appropriate places
%         RHS = [bx; by; bz];
%
%
%         pinv = 1:length(P);
%         pinv(P) = 1:length(P);
%
%         PRHS = full(RHS); PRHS = PRHS(P);
%
%         [myL, myD] = ldlsparse(L, P);
%         myL = myL + eye(size(myL));
%
%         % first solve
%         u = conelp_forwardsub(myL,PRHS);
%         v = conelp_byDiag(myD,u);
%         w = conelp_backwardsub(myL',v);
%
%         dx(P,1) = w;
%
%
%
%         % iterative refinement
%
%         for i = 1:3
%
%             % variables
%             x = dx(1:n);
%             y = dx(n+1:n+p);
%             z = dx(n+p+1:end);
%
%             % errors
%             ex = bx - A'*y - Gtilde'*z;
%             ey = by - A*x;
%             ez = bz - Gtilde*x + V*z;
%
%             % permute error
%             e(pinv,1) = [ex; ey; ez];
%
%             % solve for correction
%             u = conelp_forwardsub(myL,e);
%             v = conelp_byDiag(myD,u);
%             Pdx = conelp_backwardsub(myL',v);
%
%             if( any(isnan(Pdx)) )
%                 error('Linear solver returned NaN in ITERATIVE REFINEMENT STEP');
%             end
%             dx = dx + Pdx(pinv);
%
%         end
%
%         % copy out variables
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         z = dx(n+p+1:end);
%
%     case 'ldlsparse2'
%
%         % prepare RHS
%         %RHS = [bx; by; W\bz];
%         RHS = [bx; by; Winv*bz];
%
%         % permute
%         pinv = 1:length(P);
%         pinv(P) = 1:length(P);
%         PLP = L(P,P);
%         PRHS = full(RHS); PRHS = PRHS(P);
%         Pdx = ldlsparse (PLP, [], PRHS);
%         dx = Pdx(pinv);
%         PKP = K(P,P);
%
%         % iterative refinement
%         %         bnorm = norm(PRHS, inf);
%         for i = 1:20
%             e = PRHS - PKP*dx(P);
%             %             errnorm = norm(e, inf);
%             %             if( errnorm/bnorm < 1e-14 ), break; end
%             Pdx = ldlsparse(PLP, [], e);
%             dx = dx + Pdx(pinv);
%         end
%         %         fprintf('It.Ref.Steps: %d\n',i-1);
%
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         %z = W'\(W\G*x - W\bz);
%         %z = Winv*(Gtilde*x - Winv*bz);
%         z = Winv*dx(n+p+(1:m));
%
%     case 'ldlsparse3'
%
%         %         % prepare RHS
%         %         %RHS = [bx; by; W\bz];
%         RHS = [Winv*bz; bx; by];
%
%         % permute
%         pinv = 1:length(P);
%         pinv(P) = 1:length(P);
%         PLP = L(P,P);
%         PRHS = full(RHS); PRHS = PRHS(P);
%         Pdx = ldlsparse (PLP, [], PRHS);
%         dx = Pdx(pinv);
%         PKP = K(P,P);
%
%         % iterative refinement
%         %         bnorm = norm(PRHS, inf);
%         for i = 1:3
%             e = PRHS - PKP*dx(P);
%             %             errnorm = norm(e, inf);
%             %             if( errnorm/bnorm < 5e-15 ), break; end
%             Pdx = ldlsparse(PLP, [], e);
%             dx = dx + Pdx(pinv);
%         end
%         %         fprintf('It.Ref.Steps: %d\n',i);
%
%         % dx = L \ RHS;
%
%         x = dx(m+(1:n));
%         y = dx(m+n+1:end);
%         %z = W'\(W\G*x - W\bz);
%         %z = Winv*(Gtilde*x - Winv*bz);
%         z = Winv*dx(1:m);
%
%     case 'ldlsparse5'
%
%         % prepare RHS with zeros at appropriate places
%         bztilde = bz(1:dims.l);
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + (1:dims.q(k));
%             bztilde = [bztilde; bz(coneidx); 0; 0];
%         end
%         RHS = [bx; by; bztilde];
%
%         % permutation
%         pinv = 1:length(P);
%         pinv(P) = 1:length(P);
%         PRHS = full(RHS); PRHS = PRHS(P);
%
%         % factor
%         [myL, myD] = ldlsparse(L, P);
%         myL = myL + eye(size(myL));
%
%         % solve: forward sub, divide by D, backward sub
%         u = conelp_forwardsub(myL,PRHS);
%         v = conelp_byDiag(myD,u);
%         w = conelp_backwardsub(myL',v);
%         dx(P,1) = w;
%
%         % iterative refinement
%         for i = 1:3
%
%                 % variables
%                 x = dx(1:n);
%                 y = dx(n+1:n+p);
%                 z = dx(n+p+1:end);
%
%                 % errors
%                 ex = bx - A'*y - Gtilde'*z;
%                 ey = by - A*x;
%                 ez = bztilde - Gtilde*x + V*z;
%
%                 % permute error
%                 e(pinv,1) = [ex; ey; ez];
%
%                 % solve
%                 u = conelp_forwardsub(myL,e);
%                 v = conelp_byDiag(myD,u);
%                 Pdx = conelp_backwardsub(myL',v);
%
%                 if( any(isnan(Pdx)) )
%                     error('Linear solver returned NaN in ITERATIVE REFINEMENT STEP');
%                 end
%                 dx = dx + Pdx(pinv);
%         end
%
%         % copy out variables
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         z = dx(n+p+1:n+p+dims.l);
%         for k = 1:length(dims.q)
%             coneidx = dims.l + sum(dims.q(1:k-1)) + 2*(k-1) + (1:dims.q(k));
%             z = [z; dx(n+p+coneidx)];
%         end
%
%     case 'ldlsparse6'
%
%         if( ~exist('ppp','var') )
%             ppp = zeros(size(L,1),1);
%         end
%
%         % prepare RHS
%         RHS = [bx; by; bz];
%         PRHS = full(RHS(P));
%
%         % 1st solve
%         %[LL,DD] = ldlsparse(sparse(K-ppp*ppp'), P);
%         %LL = LL + eye(size(LL));
%         LL = L;
%         DD = D;
%         for i = 1:size(ppp,2)
%             [LL,DD] = rank1update(LL,DD,-1,ppp(P,i));
%         end
%
% %         u = conelp_forwardsub(LL,PRHS);
% %         v = conelp_byDiag(DD,u);
% %         dx(P,1) = conelp_backwardsub(LL',v);
%           u = LL\PRHS; v = DD\u; Pdx = LL'\v;
%           dx(P,1) = Pdx;
% %         dx(P,1) = conelp_forwardsub_bydiag_backwardsub_r1(L,D,-1,v(P),PRHS);
%
%         % iterative refinement
%         for i = 1:5
%
%             % variables
%             x = dx(1:n);
%             y = dx(n+1:n+p);
%             z = dx(n+p+1:end);
%
%             % errors
%             ex = bx - A'*y - G'*z; fprintf('norm(ex): %4.2e  ',norm(ex));
%             ey = by - A*x;         fprintf('norm(ey): %4.2e  ',norm(ey));
%             ez = bz - G*x + V*z;   fprintf('norm(ez): %4.2e\n',norm(ez));
%
%
%             % error
%             e = [ex; ey; ez];
%
%             % solve & refine
% %             dx(P,1) = dx(P,1) + conelp_forwardsub_bydiag_backwardsub_r1(L,D,-1,v(P),e(P));
%             u = conelp_forwardsub(LL,e(P));
%             v = conelp_byDiag(DD,u);
%             dx(P,1) = dx(P,1) + conelp_backwardsub(LL',v);
%         end
%         fprintf('--\n');
%
%         % variables
%         x = dx(1:n);
%         y = dx(n+1:n+p);
%         z = dx(n+p+1:end);
%
%     case 'conelp_ldl'
%         u = L\RHS; v = D\u; dx = L'\v;
%
%     otherwise,
%         error('Unknown linear solver %s, exiting',LINSOLVER);
% end

