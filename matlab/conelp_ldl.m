function [L,D] = conelp_ldl(A,n,m,p)
% Naive LDL factorization method. For testing purposes only.

[n1,m1] = size(A);
assert(n1==m1, 'Matrix must be square');

d = NaN(1,n1);
L = eye(n1);

% delta = 1e-7;

for j=1:n1
%     temp = 0;
%     for k=1:j-1
%         temp = temp + L(j,k)^2*d(k);
%     end
    temp = sum(L(j,1:j-1).^2.*d(1:j-1));
    d(j) = A(j,j) - temp;
    
    
    for i = j+1:n1
        temp = sum(L(i,1:j-1).*L(j,1:j-1).*d(1:j-1));
        L(i,j) = (A(i,j) - temp) / d(j);
    end
end

D = diag(d);