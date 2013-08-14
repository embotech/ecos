clear;

dirs = dir;
k = 1;

excludedir = {'docuex'};

for dd = 1:length(dirs)
    if( ~isempty(strfind(dirs(dd).name,'.')) )
        continue;
    end
    
    if( any(strcmpi(excludedir,dirs(dd).name)) )
        continue;
    end
    
%     if( isempty(strfind(dirs(dd).name,'least_squares')) )
%         continue;
%     end 
    
    fprintf('Testing ''%s'' ', dirs(dd).name);
    cd(dirs(dd).name);
    
    % test whether this test is randomized or not
    fid = fopen('data.m');
    datdef = textscan(fid,'%s');
    fclose(fid);
    N = 1;
    for kt = 1:length(datdef{1})
        temp = strfind(datdef{1}(kt),'rand');
        if( ~isempty(temp{1}) )
            fprintf('on randomized data ');
            N = 1e5;
            break;
        end
    end
    fprintf('...\n');
    
    
    delStr = '';
    for i = 1:N
        data; evalc('ecos_solver');
        exitflags(k,i) = info_.exitflag;
        
        bar = repmat(sprintf('='),1,round(i/N*30));
        msg = sprintf('[%-30s] %d/%d', bar, i, N);
        fprintf([delStr, msg]);
        delStr = repmat(sprintf('\b'), 1, length(msg));
    end
    cd ..
    fprintf('\n');
    succ(k) = sum(exitflags(k,1:N) >= 0);
    fprintf('%d (%4.2f%%) converged to a solution within desired accuracy.\n\n', succ(k), succ(k)/N*100);
    
    k = k+1;
    
end
