function dumpMatrix(x,name)

fid = fopen(name, 'w');

[m, n] = size(x);

for i = 1:m
    for j = 1:n
        fprintf(fid,'%+18.16e',x(i,j));
        if( j < n )
            fprintf(fid,'\t');
        end
    end
    fprintf(fid,'\n');
end

fclose(fid);

fprintf('Matrix of size %dx%d ''%s'' successfully written.\n', m, n, name);