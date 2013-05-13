close all; clear all;
% variable sizes
sizes = [3 10 30 100 300];% 3000 10000];
cone_sizes = [0 2 4 8 16 32 64] + 1;

svm_objval = [];
port_objval = [];
%results_svm = cell(length(sizes), length(sizes), length(cone_sizes));
%results_portfolio = cell(length(sizes), length(sizes), length(cone_sizes));
for i = 1:length(sizes),
    m = sizes(i);
    for j = 1:length(sizes),
        n = sizes(j);
        
        X = 2*randn(m,n);
        Y = -2*randn(m,n);
        gamma = 10;

        mu = 2*exp(randn(n,1));
        F = 10*randn(n,m);
        D = diag(rand(n,1));
        
        for k = 1:length(cone_sizes),
            cone_size = cone_sizes(k);
                                
            fprintf('!!!!!!!!!!!!\n');
            fprintf('! Running with m = %d, n = %d, cone size = %d\n', m,n,cone_size);
            fprintf('!!!!!!!!!!!!\n');
            fprintf('\n');
            
            if cone_size < 3
                command = sprintf('!python create_matlab_source.py -m %d -n %d --svm', m,n);
                %command_cvx = sprintf('!python create_matlab_source.py -m %d -n %d --svm --cvx', m,n);
            else
                command = sprintf('!python create_matlab_source.py -m %d -n %d -q %d --svm', m,n, cone_size);
                %command_cvx = sprintf('!python create_matlab_source.py -m %d -n %d -q %d --svm --cvx', m,n, cone_size);
            end
            f1 = evalc(command);
            

            try
            s1 = evalc(f1);
%             t = cvx_optval;
%             a1 = a;
%             b1 = b;
%             x1 = x;
            result.m = m; result.n = n;
            result.cone_size = cone_size;
            result.iter = info.iter;
            result.runtime = info.timing.runtime;
            results_svm(i,j,k) = result;
            catch e
                f1
                disp(e)
                return
            end
            svm_objval(i,j,k) = info.pcost;
            

%             f2 = evalc(command_cvx);
%             f2
%             eval(f2)
%             [a a1]
%             [b b1]
%             [cvx_optval t]
%             return

            if cone_size < 3
                command = sprintf('!python create_matlab_source.py -m %d -n %d --portfolio', m,n);
            else
                command = sprintf('!python create_matlab_source.py -m %d -n %d -q %d --portfolio', m,n, cone_size);
            end
            f2 = evalc(command);


            try
            s2 = evalc(f2);
            result.m = m; result.n = n;
            result.cone_size = cone_size;
            result.iter = info.iter;
            result.runtime = info.timing.runtime;
            results_portfolio(i,j,k) = result;
            catch e
                f2
                disp(e)
                return
            end
            port_objval(i,j,k) = info.pcost;
        end
    end
end

%%
port_fillin = [];
svm_fillin = [];
port_nnzL = [];
port_nnzK = [];
svm_nnzL = [];
svm_nnzK = [];
port_q = [];
svm_q = [];
port_iter = [];
svm_iter = [];
port_time = [];
svm_time = [];
for i = 1:length(results_portfolio(:)),
    if ~isempty(results_portfolio(i).fillin)
        port_fillin(end+1,:) = results_portfolio(i).fillin;
        port_nnzK(end+1,:) = results_portfolio(i).nnzK;
        port_nnzL(end+1,:) = results_portfolio(i).nnzL;
        port_iter(end+1) = results_portfolio(i).iter;
        port_time(end+1) = results_portfolio(i).runtime;
        port_q(end+1) = results_portfolio(i).cone_size;
    end
end


for i = 1:length(results_svm(:)),    
    if ~isempty(results_svm(i).fillin)
        svm_fillin(end+1,:) = results_svm(i).fillin;
        svm_nnzK(end+1,:) = results_svm(i).nnzK;
        svm_nnzL(end+1,:) = results_svm(i).nnzL;
        svm_iter(end+1) = results_svm(i).iter;
        svm_time(end+1) = results_svm(i).runtime;
        svm_q(end+1) = results_svm(i).cone_size;
    end
end

s = [];
e = [];
U = [];
L = [];
ennz = [];
Unnz = [];
Lnnz = [];

e_iter = [];
s_iter = [];
e_time = [];
s_time = [];
for i = 1:length(cone_sizes),
    mask = (svm_q == cone_sizes(i));
    s(end+1) = cone_sizes(i);
    e(end+1,:) = mean(svm_fillin(mask,:));
    U(end+1,:) = max(svm_fillin(mask,:));
    L(end+1,:) = min(svm_fillin(mask,:));
    
    e_iter(end+1) = mean(svm_iter(mask));
    s_iter(end+1) = std(svm_iter(mask));
    e_time(end+1) = mean(svm_time(mask));
    s_time(end+1) = std(svm_time(mask));
    
    ennz(end+1,:) = [mean(svm_nnzK(mask,:)) mean(svm_nnzL(mask,:))];
    Unnz(end+1,:) = [max(svm_nnzK(mask,:)) max(svm_nnzL(mask,:))];
    Lnnz(end+1,:) = [min(svm_nnzK(mask,:)) min(svm_nnzL(mask,:))];
end

figure(1)
h = bar(e);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('fill in');
legend('expanded', 'dense', 'Location', 'BestOutside', 'Orientation', 'horizontal');
hold on
title('SVM Problem');

numgroups = size(e,1);
numbars = size(e,2);

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, e(:,i), e(:,i) - L(:,i), U(:,i) - e(:,i), 'k', 'linestyle', 'none');
end
hold off

print -depsc svm_fillin.eps

figure(2)
h = bar(ennz);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('nnz');
legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

print -depsc svm_nnz.eps

figure(3)
h = stem(e_iter);
set(gca, 'XTickLabel', s);
ylabel('iters');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

print -depsc svm_iter.eps

figure(4)
h = stem(e_time);
set(gca, 'XTickLabel', s);
ylabel('time (s)');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');


print -depsc svm_time.eps


s = [];
e = [];
U = [];
L = [];
ennz = [];
Unnz = [];
Lnnz = [];

e_iter = [];
s_iter = [];
e_time = [];
s_time = [];

for i = 1:length(cone_sizes),
    mask = (port_q == cone_sizes(i));
    s(end+1) = cone_sizes(i);
    e(end+1,:) = mean(port_fillin(mask,:));
    U(end+1,:) = max(port_fillin(mask,:));
    L(end+1,:) = min(port_fillin(mask,:));
    
    e_iter(end+1) = mean(port_iter(mask));
    s_iter(end+1) = std(port_iter(mask));
    e_time(end+1) = mean(port_time(mask));
    s_time(end+1) = std(port_time(mask));
    
    ennz(end+1,:) = [mean(port_nnzK(mask,:)) mean(port_nnzL(mask,:))];
    Unnz(end+1,:) = [max(port_nnzK(mask,:)) max(port_nnzL(mask,:))];
    Lnnz(end+1,:) = [min(port_nnzK(mask,:)) min(port_nnzL(mask,:))];
end

figure(5)
h = bar(e);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('fill in');
legend('expanded', 'dense', 'Location', 'BestOutside', 'Orientation', 'horizontal');
hold on
title('Portfolio Problem');

numgroups = size(e,1);
numbars = size(e,2);

groupwidth = min(0.8, numbars/(numbars+1.5));

for i = 1:numbars
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
errorbar(x, e(:,i), e(:,i) - L(:,i), U(:,i) - e(:,i), 'k', 'linestyle', 'none');
end
hold off

print -depsc port_fillin.eps


figure(6)
h = bar(ennz);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('nnz');
legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

print -depsc port_nnz.eps

figure(7)
h = stem(e_iter);
set(gca, 'XTickLabel', s);
ylabel('iters');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');


print -depsc port_iter.eps

figure(8)
h = stem(e_time);
set(gca, 'XTickLabel', s);
ylabel('time (s)');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

print -depsc port_time.eps

