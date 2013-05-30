close all; clear all
tests = dir('DIMACS/*.mat');

% optvals taken from
% http://dimacs.rutgers.edu/Challenges/Seventh/Instances/tablestat.html
objval = struct(...
    'nql30', -0.9460, ...
    'nql60', -0.935, ...
    'nql180', 0, ... % n/a
    'qssp30', -6.4966749, ...
    'qssp60', -6.5627049, ...
    'qssp180', 0, ... % n/a
    'nb', -0.05070309, ...
    'nb_L1', -13.012337, ...
    'nb_L2', -1.62897198, ...
    'nb_L2_bessel', -0.102569511, ...
    'sched_50_50_orig', 26673.00, ...
    'sched_100_50_orig', 181889.9, ...
    'sched_100_100_orig', 717367, ...
    'sched_200_100_orig', 141360.4464, ...
    'sched_50_50_scaled',7.8520384, ...
    'sched_100_50_scaled', 6.716503, ...
    'sched_100_100_scaled', 27.3307, ...
    'sched_200_100_scaled', 51.81196099 ...
);


for i = 1:length(tests),
    clear A At b c K
    test_name = tests(i).name;
    test_name = test_name(1:end-4); % strip .mat
    f = ['DIMACS/' tests(i).name];
    disp(['solving ' f]);
    load(f)
    
    [m1,n1] = size(b);
    [m2,n2] = size(c);
    if m1 == 1,
        b = b';
    end
    if m2 == 1,
        c = c';
    end
    if (exist('A','var'))
        data = struct('A', sparse(A'), 'b', full(c), 'c', -full(b));
    else
        data = struct('A', sparse(At), 'b', full(c), 'c', -full(b));
    end
    cone = struct('f', 0, 'q', K.q', 'l', K.l);
    
    [x_m, y_m, info] = ecos(data.c, data.A, data.b, cone);
    
    err.(test_name) = abs( data.c'*x_m + objval.(test_name) ) ;
    test_info.(test_name) = info;

    
%     [m,n] = size(data.A);
%     
%     ns_soc = cone.q;
%     idxs = [cone.f+cone.l; cone.f+cone.l + cumsum(ns_soc)];
%     
%     cvx_begin
%         cvx_solver sedumi
%         variable xcvx(n)
%         variable scvx(m)
% 
% 
%         minimize (data.c'*xcvx)
%             data.A*xcvx + scvx == data.b
%             scvx(1:cone.f) == 0
%             scvx(cone.f+1:cone.f+cone.l) >= 0
%             for kk =1:cone.k_soc
%                 norm(scvx(idxs(kk) + 2: idxs(kk) + ns_soc(kk))) <= scvx(idxs(kk)+1)
%             end
%     cvx_end
end