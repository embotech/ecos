function result = fillintest(A,G,dims)

oldpath = addpath('./tools', './tools/LDL/MATLAB');

[p,n] = size(A);

eps = 1e-14;
delta = 7e-14;

spyplot = 1;

%% 1. what ecos does
Gtilde = conelp_stretch(G,dims,2);
mtilde = size(Gtilde,1);
nvars = n+p+mtilde;

% use random scaling
V = eye(dims.l);
W = eye(dims.l);
for i = 1:length(dims.q)
    [Vi,Wi] = randomScaling(dims.q(i));
    V = blkdiag(V,sparse(Vi));
    W = blkdiag(W,sparse(Wi));
    i
end

% build KKT matrix
K = conelp_KKTmatrix(A,Gtilde,V,0);
nnzK = nnz(u(K));

% get permutation
P = conelp_getPerm(K~=0);

% get signs
S = [ones(n,1); -ones(p,1); -ones(dims.l,1)];
for k = 1:length(dims.q)
    S = [S; -ones(dims.q(k),1); -1; 1];
end

% factor
[L,D] = sldlsparse(sparse(K),P,S,eps,delta);
assert( all(all(~isnan(L))),'L contains NaNs' );
nnzL = nnz(L);% + nvars; 

% report
fprintf('Fill-in report for ECOS with sparse scalings:\n');
fprintf('\nNumber of non-zeros in K: %d\n', nnzK);
fprintf('Number of non-zeros in L: %d\n', nnzL);
fillin(1)=nnzL/nnzK;
fprintf('<==> fill-in fatcor: *%4.2f*\n\n', fillin(1));

% plot
if( spyplot )
    figure;
    subplot(2,2,1); spy(K);
    title('KKT (L sparse SOC blocks)');
    subplot(2,2,2); spy(L);
    title('L (sparse SOC blocks)');
end

nnzs(1) = nnzK;

%% this would happen if we used the non-expanded version
m = size(G,1);
nvars = n+p+m;

% KKT matrix
K = conelp_KKTmatrix(A,G,W,0);
nnzK = nnz(u(K));

% permutation
P = conelp_getPerm(K~=0);

% signs
S = [ones(n,1); -ones(p,1); -ones(m,1)];

% factor
[L,D] = sldlsparse(sparse(K),P,S,eps,delta);
assert( all(all(~isnan(L))),'L contains NaNs' );
nnzL = nnz(L);% + nvars;

% report
fprintf('Fill-in report for ECOS with dense scalings:\n');
fprintf('\nNumber of non-zeros in K: %d\n', nnzK);
fprintf('Number of non-zeros in L: %d\n', nnzL);
fillin(2)=nnzL/nnzK;
fprintf('<==> fill-in factor: *%4.2f*\n\n', fillin(2));

% plot
if( spyplot )
    subplot(2,2,3);
    spy(K);
    title('KKT (dense SOC blocks)');
    subplot(2,2,4);
    spy(L);
    title('L (dense SOC blocks)'); 
end

nnzs(2) = nnzK;

result.fillin = fillin;
result.nnz = nnzs;

addpath(oldpath)    % restore