function txt = cg_dump_conelpproblem(fn,c,G,h,dims,A,b,prefix)
% Dump a conic problem in sparse format to a text file - for use with ECOS.
% USAGE:
%
%   txt = cg_dump_conelpproblem(c,G,h,dims) dumps the problem in sparse
%   format such that it can be directly used with the C interface of ECOS.
%   The output argument TXT is the text that has been created in C
%   declaration format.
%
%   txt = cg_dump_conelpproblem(c,G,h,dims,prefix) same as above, but
%   prepends variable names with the string prefix.
%
%
% If in addition equality constraints are present you can use the following
% commands:
%
%   txt = cg_dump_conelpproblem(c,G,h,dims,A,b) same as above, with
%   additional problem parameters A and b.
%
%   txt = cg_dump_conelpproblem(c,G,h,dims,A,b,prefix) to prepend variables
%   with the string prefix.
%
%
% To directly save the dumped data into a file, use:
%
%   txt = cg_dump_conelpproblem(fn,...) additionally saves TXT to the file fn.
%
% (c) Alexander Domahidi, embotech GmbH, March 2014. <domahidi@embotech.com>


%% check which version we're running
% we need at least 4 arguments: c,G,h,dims
if nargin < 4
    error('At least 4 arguments needed: c, G, h, dims');
end

% user supplied c, G, h, dims in form of fn, c, G, h
if nargin == 4
    assert( ~ischar(fn), 'At least 5 arguments needed: fn, c, G, h, dims' );
    dims = h;
    h = G;
    G = c;
    c = fn;
    A = [];
    b = [];
    prefix = '';
    savefile = false;
end

% user supplied either (c, G, h, dims, prefix) or (fn, c, G, h, dims)
if( nargin == 5 )
    if( ischar(fn) )
        % user supplied (fn, c, G, h, dims)
        prefix = '';
        savefile = true;
    else
        % user supplied (c, G, h, dims, prefix) as (fn, c, G, h, dims)
        % have to play rochade
        assert(ischar(dims),'Invalid arguments: c, G, h, dims, prefix expected - 5th argument must be a string');
        prefix = dims;
        dims = h;
        h = G;
        G = c;
        c = fn;
        savefile = false;
    end
    A = [];
    b = [];
end


% user supplied (c,G,h,dims,A,b) as (fn,c,G,h,dims,A)
if( nargin == 6 )
    assert( ~ischar(fn),'6 arguments supplied - c,G,h,dims,A,b expected but c is string' );
        % move fn,c,G,h,   dims,A to 
        %       c,G,h,dims,  A ,b
        b = A;
        A = dims;
        dims = h;
        h = G;
        G = c;
        c = fn;
        prefix = '';
end

% user supplied either (c,G,h,dims,A,b,prefix) as (fn,c,G,h,dims,A,b)
%                   or (fn,c,G,h,dims,A,b) as (fn,c,G,h,dims,A,b)
if( nargin == 7 )
    if( ischar(fn) )
        % user supplied (fn,c,G,h,dims,A,b) as (fn,c,G,h,dims,A,b)
        prefix = '';
        savefile = true;
    else
        % user supplied (c,G,h,dims,A,b,prefix) as (fn,c,G,h,dims,A,b)
        assert( ischar(b), '7 arguments provided - expected c,G,h,dims,A,b,prefix but prefix is not a string');
        prefix = b;
        b = A;
        A = dims;
        dims = h;
        h = G;
        G = c;
        c = fn;
        savefile = false;
    end
end

%% do the actual codegen
[m,n] = size(G);
p = size(A,1);
txt = {};
txt = [txt; 'idxint ', prefix, 'n = ',num2str(n),';'];
txt = [txt; 'idxint ', prefix, 'm = ',num2str(m),';'];
txt = [txt; 'idxint ', prefix, 'p = ',num2str(p),';'];
txt = [txt; 'idxint ', prefix, 'l = ',num2str(dims.l),';'];
txt = [txt; 'idxint ', prefix, 'ncones = ',num2str(length(dims.q)),';'];

disp('Creating string for c...');
txt = [txt; cg_mat2c('pfloat',[prefix,'c'],c,'%20.18e')];

disp('Creating string for G...');
txt = [txt; cg_dump_spmat(G,[prefix,'G'],'pfloat','idxint')];

disp('Creating string for h...');
txt = [txt; cg_mat2c('pfloat',[prefix,'h'],h,'%20.18e')];

disp('Creating string for dims...');
txt = [txt; cg_mat2c('idxint',[prefix,'q'],dims.q,'%d')];

disp('Creating string for A...');
txt = [txt; cg_dump_spmat(A,[prefix,'A'],'pfloat','idxint')];

disp('Creating string for b...');
txt = [txt; cg_mat2c('pfloat',[prefix,'b'],b,'%20.18e')];
    
% write file
if( savefile )    
    disp('Flushing string to file...');
    cg_dumpfile(fn,txt);
end