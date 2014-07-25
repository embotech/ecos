#A = [0 1; 1 0]

[ir, row, pr] = find(A);
ir = ir - 1;
[~, jc] = unique(row);
jc = [0; jc];

jc = jc'
ir = ir'
pr = pr'


