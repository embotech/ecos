function txt = cg_dump_conelpproblem(fn,c,G,h,dims,A,b)
% Dump a conic problem in sparse format to a text file.
% USAGE:
%
%   txt = cg_dump_conelpproblem(fn,c,G,h,dims)
%
% or
%
%   txt = cg_dump_conelpproblem(fn,c,G,h,dims,A,b)
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.


if( nargin == 5 )
    A = [];
    b = [];
end
assert(ischar(fn),'First argument must be a string containing the file name');

[m,n] = size(G);
p = size(A,1);
txt = {};
txt = [txt; 'idxint n = ',num2str(n),';'];
txt = [txt; 'idxint m = ',num2str(m),';'];
txt = [txt; 'idxint p = ',num2str(p),';'];
txt = [txt; 'idxint l = ',num2str(dims.l),';'];
txt = [txt; 'idxint ncones = ',num2str(length(dims.q)),';'];

disp('Creating string for c...');
txt = [txt; cg_mat2c('pfloat','c',c,'%20.18e')];

disp('Creating string for G...');
txt = [txt; cg_dump_spmat(G,'G','pfloat','idxint')];

disp('Creating string for h...');
txt = [txt; cg_mat2c('pfloat','h',h,'%20.18e')];

disp('Creating string for dims...');
txt = [txt; cg_mat2c('idxint','q',dims.q,'%d')];

disp('Creating string for A...');
txt = [txt; cg_dump_spmat(A,'A','pfloat','idxint')];

disp('Creating string for b...');
txt = [txt; cg_mat2c('pfloat','b',b,'%20.18e')];
    

% write file
disp('Flushing string to file...');
cg_dumpfile(fn,txt);