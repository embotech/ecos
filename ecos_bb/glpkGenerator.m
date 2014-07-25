#A = [ 2  1
#      3  4 
#      1  0
#     -1  0
#      0 -1];

#b = [ 4; 12; 1; -1; 0];
#c = [-1, -1];

c = [-3 -5 6 9 10 -10];
A = [-2 6 -3 4 1 -2;
     -5 -3 1 3 -2 1;
     5 -1 4 -2 2 -1];
A = -A;
b = -[2;-1;3];

ctype = 'UUU';
vartype = 'IIIIII';

sparse(A)
lb = zeros(1,6);
ub = ones(1,6);

[xopt, fmin] = glpk (c, A, b, lb, ub, ctype, vartype, 1)