close all; clear all;
% variable sizes
%sizes = [3 10 30 100 300];% 3000 10000];
svm_size = [1000 30];
port_size = [30 1000];
N = 1;
cone_sizes = [0 2 4 8 16 32 64] + 1;

svm_objval = [];
port_objval = [];
%results_svm = cell(length(sizes), length(sizes), length(cone_sizes));
%results_portfolio = cell(length(sizes), length(sizes),
%length(cone_sizes));

m = svm_size(1); n = svm_size(2);

X = 2*randn(m,n);
Y = -2*randn(m,n);
gamma = 10;

m = port_size(1); n = port_size(2);
mu = 2*exp(randn(n,1));
F = 10*randn(n,m);
D = diag(rand(n,1));

for k = 1:length(cone_sizes),
    cone_size = cone_sizes(k);

    m = svm_size(1); n = svm_size(2);

    fprintf('!!!!!!!!!!!!\n');
    fprintf('! Running SVM with m = %d, n = %d, cone size = %d\n', m,n,cone_size);
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
    results_svm(k) = result;
    catch e
        f1
        disp(e)
        return
    end
    svm_objval(k) = info.pcost;


%             f2 = evalc(command_cvx);
%             f2
%             eval(f2)
%             [a a1]
%             [b b1]
%             [cvx_optval t]
%             return

    m = port_size(1); n = port_size(2);

    fprintf('!!!!!!!!!!!!\n');
    fprintf('! Running PORT with m = %d, n = %d, cone size = %d\n', m,n,cone_size);
    fprintf('!!!!!!!!!!!!\n');
    fprintf('\n');

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
    results_portfolio(k) = result;
    catch e
        f2
        disp(e)
        return
    end
    port_objval(k) = info.pcost;
end

%%

print_figs = 0;

s = [];
e = [];
ennz = [];

e_iter = [];
e_time = [];
for i = 1:length(cone_sizes),
    s(end+1) = cone_sizes(i);
    e(end+1,:) = results_svm(i).fillin;
    
    e_iter(end+1) = results_svm(i).iter;
    e_time(end+1) = results_svm(i).runtime;
    
    ennz(end+1,:) = [results_svm(i).nnzK results_svm(i).nnzL];
end

figure(1)
h = bar(e);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('fill in');
legend('expanded', 'dense', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

if print_figs
print -depsc svm_fillin.eps
end

figure(2)
h = bar(ennz);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('nnz');
legend('K', 'K', 'K', 'dense K', 'AMD', 'SYMAMD', 'METIS', 'dense METIS', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

if print_figs
print -depsc svm_nnz.eps
end

figure(3)
h = stem(e_iter);
set(gca, 'XTickLabel', s);
ylabel('iters');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

if print_figs
print -depsc svm_iter.eps
end

figure(4)
h = stem(e_time);
set(gca, 'XTickLabel', s);
ylabel('time (s)');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('SVM Problem');

if print_figs
print -depsc svm_time.eps
end

s = [];
e = [];
ennz = [];

e_iter = [];
e_time = [];
for i = 1:length(cone_sizes),
    s(end+1) = cone_sizes(i);
    e(end+1,:) = results_portfolio(i).fillin;
    
    e_iter(end+1) = results_portfolio(i).iter;
    e_time(end+1) = results_portfolio(i).runtime;
    
    ennz(end+1,:) = [results_portfolio(i).nnzK results_portfolio(i).nnzL];
end

figure(5)
h = bar(e);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('fill in');
legend('expanded', 'dense', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

if print_figs
print -depsc port_fillin.eps
end

figure(6)
h = bar(ennz);
set(h, 'BarWidth', 1);
set(gca, 'XTickLabel', s);
ylabel('nnz');
legend('K', 'K', 'K', 'dense K', 'AMD', 'SYMAMD', 'METIS', 'dense METIS', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

if print_figs
print -depsc port_nnz.eps
end

figure(7)
h = stem(e_iter);
set(gca, 'XTickLabel', s);
ylabel('iters');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

if print_figs
print -depsc port_iter.eps
end

figure(8)
h = stem(e_time);
set(gca, 'XTickLabel', s);
ylabel('time (s)');
%legend('expanded KKT', 'dense KKT', 'expanded L', 'dense L', 'Location', 'BestOutside', 'Orientation', 'horizontal');
title('Portfolio Problem');

if print_figs
print -depsc port_time.eps
end
