function [L,D] = conelp_factor(K,P)
% LDL Factorization routine for CONELP solver.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% make K sparse
if( ~issparse(K) ), K = sparse(K); end

% factor & permute  
[L,D] = ldlsparse(K, P);
L = L + eye(size(L));

%% DEPRECATED 
% 
% switch( lower(LINSOLVER) )
%     
%     % MATLAB's backslash LDL factorizer
%     case 'backslash'
%         [L,D,P] = ldl(K);
%         
%     case 'backslash2'
%         [L,D,P] = ldl(K);    
%     
%     case 'backslash3'
%         [L,D,P] = ldl(K);    
%             
%     case 'ldlsparse'
%         
%         
%         
%     case 'ldlsparse4'
%         % regularize
%         eps = 1e-7;      
%         K(1:n,1:n) = eps*speye(n);
%         K(n+(1:p),n+(1:p)) = -eps*speye(p);
%  
%         % factor - we don't do anything here yet, factorization is done
%         % with the solve for the moment
%         L = K;
%         D = [];
%         P = amd(K);
%         
%     case 'ldlsparse5'
%         % regularize
%         eps = 1e-7;      
%         K(1:n,1:n) = eps*speye(n);
%         K(n+(1:p),n+(1:p)) = -eps*speye(p);
%  
%         % factor - we don't do anything here yet, factorization is done
%         % with the solve for the moment
%         L = K;
%         D = [];
%         P = amd(K);
%        
%     case 'ldlsparse2'
%         % regularize
%         eps = 1e-7;        
%         K(1:n,1:n) = +eps*eye(n);
%         K(n+(1:p),n+(1:p)) = -eps*eye(p);
%         K = sparse(K);
%         L = K;
%         D = [];
%         P = amd(K);   
%         
%     case 'ldlsparse3'
%         L = K;
% %         regularize
%         eps = 1e-8;        
%         L(m+(1:n),m+(1:n)) = +eps*eye(n);
%         L(m+n+(1:p),m+n+(1:p)) = -eps*eye(p);
%         L = sparse(L);
%         D = [];
%         P = amd(L);   
%         
%     case 'ldlsparse6'
%         % factor
%         if( ~exist('P','var'))
%             P = amd(K);
%         end
%         PKP = sparse(K(P,P));
%         [L,D] = ldlsparse(PKP);
%         L = L + speye(size(L));          
%    
%         
%     case 'conelp_ldl'
%         %         eps = 1e-5;
%         %         K(1:n,1:n) = +eps*speye(n);
%         %         K(n+(1:p),n+(1:p)) = -eps*speye(p);
%         [L,D] = conelp_ldl(K,n,m,p);
%         P = eye(size(K));
%         
%     otherwise,
%         error('Unknown linear solver %s, exiting',LINSOLVER);
% end
