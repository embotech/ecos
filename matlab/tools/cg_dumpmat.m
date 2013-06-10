function matstr = cg_dumpmat(M,elsep,fstr)
% CG_DUMPMAT Convert matrix to string.
%    MATSTR = CG_DUMPMAT(M) converts the matrix M to a string which is
%    returned in MATRSTR. The elements of M will be seperated by commas.
%
%    MATSTR = CG_DUMPMAT(M,ELSEP) converts the matrix M to a string which
%    is returned in MATRSTR. The elements of M will be seperated by the
%    string ELSEP.
%
%    MATSTR = CG_DUMPMAT(M,ELSEP,FORMATSTR) converts the matrix M to a
%    string with the elements of M seperated by string ELSEP and formatted
%    according to FORMATSTR.
% 
% see also CG_DUMPFILE

dimi = size(M,1);
dimj = size(M,2);
matstr = cell(dimi,1);

if( nargin < 3 ) 
    fstr = '%20.18e';
end;

if( ~exist('elsep','var') )
    elsep = ',';
end

for i = 1:dimi
    if( M(i,1) == 0 )
        temp = '0.0';
    else
        temp = num2str(M(i,1),fstr);    
    end
    for j = 2:dimj
        if( M(i,j) == 0 )
            valstr = '0.0';
        else
            valstr = num2str(M(i,j),fstr);
        end
        temp = [temp,elsep,' ',valstr];
    end    
    if( i < dimi )
        temp = [temp,elsep];
    end
    matstr{i} = temp;
end

if( dimi == 1 )
    matstr = matstr{1};
end