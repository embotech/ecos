function x = conelp_byDiag(D,y)
% returns x = D \ y where D is a diagonal matrix

n = size(D,1);
x = NaN(n,1);
for i = 1:n
    assert( D(i,i) ~= 0, 'division by zero');
    x(i) = y(i) / D(i,i);
end