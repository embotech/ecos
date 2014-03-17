function txt = cg_dump_conelpproblem(fn,c,G,h,dims,A,b)
% Dump a conic problem in sparse format to a text file - for use with ECOS.
% USAGE:
%
%   txt = cg_dump_conelpproblem(c,G,h,dims) dumps the problem in sparse
%   format such that it can be directly used with the C interface of ECOS.
%   The output argument TXT is the text that has been created in C
%   declaration format.
%
%
%   txt = cg_dump_conelpproblem(c,G,h,dims,A,b) same as above, with
%   additional problem parameters A and b.
%
%
%   txt = cg_dump_conelpproblem(fn,...) additionally saves TXT to the file fn.
%
% (c) Alexander Domahidi, embotech GmbH, March 2014. <domahidi@embotech.com>


if( nargin == 4 || nargin == 5 )
    A = [];
    b = [];
end

if( nargin == 4 || nargin == 6 )
    savefile = false;
else
    savefile = true;
end
    
if( nargin == 4 )
    % move fn,c,G,h     to 
    %       c,G,h,dims
    dims = h;
    h = G;
    G = c;
    c = fn;
end

if( nargin == 6 )
    % move fn,c,G,h,   dims,A to 
    %       c,G,h,dims,  A ,b
    b = A;
    A = dims;
    dims = h;
    h = G;
    G = c;
    c = fn;
end
    
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
if( savefile )    
    disp('Flushing string to file...');
    cg_dumpfile(fn,txt);
end