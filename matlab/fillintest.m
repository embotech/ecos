function fillin = fillintest(A,G,dims)

[p,n] = size(A);


%% 1. what ecos does
Gtilde = conelp_stretch(G,dims,2);
mtilde = size(Gtilde,1);
nvars = n+p+mtilde;
Vpattern = conelp_scaling(dims, 'ldlsparse', 'pattern');
Kpattern = conelp_KKTmatrix(A,Gtilde,Vpattern,0);
nnzK = nnz(u(Kpattern));
P = conelp_getPerm(Kpattern~=0);
S = [ones(n,1); -ones(p,1); -ones(dims.l,1)];
for k = 1:length(dims.q)
    S = [S; -ones(dims.q(k),1); -1; 1];
end
[L,D] = sldlsparse(sparse(Kpattern),P,S,1e-14,7e-7);
assert( all(all(~isnan(L))),'L contains NaNs' );
nnzL = nnz(L) + nvars; 
fprintf('Fill-in report for ECOS with sparse scalings:\n');
fprintf('\nNumber of non-zeros in K: %d\n', nnzK);
fprintf('Number of non-zeros in L: %d\n', nnzL);
fillin(1)=nnzL/nnzK;
fprintf('<==> fill-in fatcor: *%4.2f*\n\n', fillin(1));

%% this would happen if we used the non-expanded version
m = size(G,1);
nvars = n+p+m;
Vpattern = conelp_scaling(dims, 'backslash', 'pattern');
Kpattern = conelp_KKTmatrix(A,G,Vpattern,0);
nnzK = nnz(u(Kpattern));
P = conelp_getPerm(Kpattern~=0);
S = [ones(n,1); -ones(p,1); -ones(m,1)];
[L,D] = sldlsparse(sparse(Kpattern),P,S,1e-14,7e-7);
assert( all(all(~isnan(L))),'L contains NaNs' );
nnzL = nnz(L) + nvars;
fprintf('Fill-in report for ECOS with dense scalings:\n');
fprintf('\nNumber of non-zeros in K: %d\n', nnzK);
fprintf('Number of non-zeros in L: %d\n', nnzL);
fillin(2)=nnzL/nnzK;
fprintf('<==> fill-in fatcor: *%4.2f*\n\n', fillin(2));