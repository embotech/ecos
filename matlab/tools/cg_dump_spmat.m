function txt = cg_dump_spmat(G,name,ftype,itype,fn)
% Dumps a sparse matrix in column compressed format into a text file.
%
%   TXT = CG_DUMP_SPMAT(G,NAME,FTYPE,ITYTPE) writes the matrix G in column  
%   compressed format in TXT such that it can be directly used by the sparse
%   matrix data structure spmat in ecos, for example. The name of the
%   variables are NAME, they are of FTYPE and the indices are ITYTPE.
%
%   Default values for FTYPE is "double" and for ITYTPE "int". If G is
%   empty, all arrays are empty and NULLed.
%
%   
%   TXT = CG_DUMP_SPMAT(G,NAME,FTYPE,ITYPE,FN) writes in addition to the 
%   above the result TXT into a file FN.
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.

if( ~exist('ftype','var') || isempty(ftype) )
    ftype = 'double';
end

if( ~exist('itype','var') || isempty(ftype) )
    ftype = 'int';
end

if( ~isempty( G ) ) 

[m,n] = size(G);
G = sparse(G);

% find jc, ir and pr
jc = NaN(1,n+1);
jc(1) = 1;
k = 1;
for j = 1:n    
    jc(j) = k;
    for i = 1:m
        v = full(G(i,j));
        if( v ~= 0 )
            ir(k) = i; 
            pr(k) = v;
            k = k+1; 
        end
    end
end
jc(n+1) = k;


% create text
txt{1,1}    = [itype,' ',name,'jc[',num2str(length(jc)),'] = {',cg_dumpmat(jc-1,',','%d'),'};'];
txt{end+1,1}= [itype,' ',name,'ir[',num2str(length(ir)),'] = {',cg_dumpmat(ir-1,',','%d'),'};'];
txt{end+1,1}= [ftype,' ',name,'pr[',num2str(length(pr)),'] = {',cg_dumpmat(pr),'};'];
% txt{end+1,1}= [itype,' ',name,'m = ',num2str(m),';'];
% txt{end+1,1}= [itype,' ',name,'n = ',num2str(n),';'];
% txt{end+1,1}= [itype,' ',name,'nnz = ',num2str(k),';'];

else
    
    % NULL all arrays since G is empty
    txt{1,1}    = [itype,' *',name,'jc = NULL;'];
    txt{end+1,1}= [itype,' *',name,'ir = NULL;'];
    txt{end+1,1}= [ftype,' *',name,'pr = NULL;'];
    
end

% dump to c file
if( nargin == 5 )
    cg_dumpfile(fn,txt);
end