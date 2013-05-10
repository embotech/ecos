close all; clear all;
% variable sizes
sizes = [3 10 30];% 100 300];% 1000 3000 10000];
cone_sizes = [0 2 4 8 16 32 64] + 1;

%results_svm = cell(length(sizes), length(sizes), length(cone_sizes));
%results_portfolio = cell(length(sizes), length(sizes), length(cone_sizes));
for i = 1:length(sizes),
    m = sizes(i);
    for j = 1:length(sizes),
        n = sizes(j);
        if m == n
            continue
        end
        for k = 1:length(cone_sizes),
            cone_size = cone_sizes(k);
        
            fprintf('!!!!!!!!!!!!\n');
            fprintf('! Running with m = %d, n = %d, cone size = %d\n', m,n,cone_size);
            fprintf('!!!!!!!!!!!!\n');
            fprintf('\n');
            
            if cone_size == 0
                command = sprintf('!python create_matlab_source.py -m %d -n %d --svm', m,n);
            else
                command = sprintf('!python create_matlab_source.py -m %d -n %d -q %d --svm', m,n, cone_size);
            end
            f1 = evalc(command);
            X = 2*randn(m,n);
            Y = -2*randn(m,n);
            gamma = 10;
            s1 = evalc(f1);
            result.m = m; result.n = n;
            result.cone_size = cone_size;
            results_svm(i,j,k) = result;



            if cone_size == 0
                command = sprintf('!python create_matlab_source.py -m %d -n %d --portfolio', m,n);
            else
                command = sprintf('!python create_matlab_source.py -m %d -n %d -q %d --portfolio', m,n, cone_size);
            end
            f2 = evalc(command);
            mu = 2*exp(randn(n,1));
            F = 10*randn(n,m);
            D = diag(rand(n,1));
            s2 = evalc(f2);
            result.m = m; result.n = n;
            result.cone_size = cone_size;
            results_portfolio(i,j,k) = result;
        end
    end
end
