close all; clear;
tests = dir('DIMACS/*.mat');

% optvals taken from
% http://dimacs.rutgers.edu/Challenges/Seventh/Instances/tablestat.html
% http://people.orie.cornell.edu/miketodd/tttmpbc.pdf
objval = struct(...
    'nql30', -0.9460, ...
    'nql60', -0.935423, ...
    'nql180', 0, ... % n/a
    'qssp30', -6.496675, ...
    'qssp60', -6.562696, ...
    'qssp180', 0, ... % n/a
    'nb', -0.05070309, ...
    'nb_L1', -13.01227, ...
    'nb_L2', -1.628972, ...
    'nb_L2_bessel', -0.102571, ...
    'sched_50_50_orig', 26673.00, ...
    'sched_100_50_orig', 181889.9, ...
    'sched_100_100_orig', 717367, ...
    'sched_200_100_orig', 141360.4464, ...
    'sched_50_50_scaled',7.852038, ...
    'sched_100_50_scaled', 67.166281, ...   % their table is wrong or something
    'sched_100_100_scaled', 27.331457, ...
    'sched_200_100_scaled', 51.812471 ...
);

k=1;
testnr = 1:length(tests)
for i = testnr
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
    
    if( ~exist('At','var') )
        At = A';
    end
    
    if( ~isfield(K,'f') ), K.f = 0; end
    
    % map sedumi primal --> ecos dual (which we give as primal)
    c_ecos = full(-b);
    G_ecos = At;
    h_ecos = full(c);
    dims = K;
    
    % solve dual with ecos
    [y_d, s_d, info_dd, z_d, x_d] = ecos(c_ecos, G_ecos, h_ecos, dims);
    ecos_fval_d = -info_dd.pcost;
    relerr_d = abs((ecos_fval_d - objval.(test_name))/objval.(test_name));
    err_d.(test_name) = relerr_d;
    iter_d.(test_name) = info_dd.iter;
    time_d.(test_name) = info_dd.timing.runtime;
    info_d.(test_name) = info_dd;
    exitflag_d.(test_name) = info_dd.exitflag;
    ecos_x_d.(test_name) = x_d;
    ecos_y_d.(test_name) = y_d;
    ecos_z_d.(test_name) = z_d;
    ecos_s_d.(test_name) = s_d;
    derr = dimacs_errors(At',b,c,K,x_d,y_d,z_d);
    pres_d.(test_name) = derr(1);
    dres_d.(test_name) = derr(2);
    if( info_dd.exitflag == 0 ), success_d(k) = 1; else success_d(k) = 0; end
    if( info_dd.exitflag == 0 )
        statstr_d.(test_name) = 'OK';
    else
        statstr_d.(test_name) = 'numerr';
    end
    fprintf('--------------------------------------------------------------------------------------------------\n');
    fprintf('DIMACS TEST ''%s'': ECOS DUAL fval: %f, true fval: %f, rel. error: %e\n', test_name, ecos_fval_d, objval.(test_name), relerr_d);
    fprintf('--------------------------------------------------------------------------------------------------\n\n\n');
    
    
    
    % solve with sedumi
    pars.fid = 0; % switch off printing
    pars.eps = 1e-6;
    pars.bigeps = 1e-6/10;
    pars.errors = 1;
    if( 1 )
        [x_s,y_s,info_ss] = sedumi(At',b,c,K,pars);
        sedumi_fval = c'*x_s;
        relerr_s = abs((sedumi_fval - objval.(test_name))/objval.(test_name));
        err_s.(test_name) = relerr_s;
        iter_s.(test_name) = info_ss.iter;
        time_s.(test_name) = sum(info_ss.timing);
        info_s.(test_name) = info_ss;
        exitflag_s.(test_name) = info_ss.numerr;
        sedumi_x_s.(test_name) = x_s;
        sedumi_y_s.(test_name) = y_s;
        derr = dimacs_errors(At',b,c,K,x_s,y_s);
        pres_s.(test_name) = derr(1);
        dres_s.(test_name) = derr(2);
%         for j = 1:length(K.q)
%             coneidx = K.f + K.l + (1:K.q(j));
%             dres_s.(test_name) = max([dres_s.(test_name), norm(ineqres(coneidx(1)) - norm(ineqres(coneidx(2:end)),2),inf)]);
%         end
        if( info_ss.numerr == 0 )
            statstr_s.(test_name) = 'OK';
        else
            statstr_s.(test_name) = 'numerr';
        end
        
    end
    
    numvar.(test_name) = size(At,1);
    numRows.(test_name) = size(At,2);
    
    k=k+1;
    
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


%% print summary in tex format
disp('DIMACS SUMMARY:');
% fprintf('Primal formulation solved %d / %d problems\n', sum(success_p), length(success_p));
fprintf('  Dual formulation solved %d / %d problems\n\n', sum(success_d), length(success_d));

for i = testnr
    test_name = tests(i).name;
    test_name = test_name(1:end-4);
    fprintf('%d & %s & %d & %d & %s & %s & %4.2e & %4.2e & %4.2e & %4.2e & %3d & %3d & %4.2f & %4.2f \\\\\n',...
        i, strrep(test_name,'_','\_'), numvar.(test_name), numRows.(test_name), ...
        statstr_s.(test_name), statstr_d.(test_name), ...
        err_s.(test_name), err_d.(test_name), ...
        pres_s.(test_name), pres_d.(test_name), ...
        iter_s.(test_name), iter_d.(test_name),...
        time_s.(test_name), time_d.(test_name) );
end


%% evalaute timings
k=0;
totalt_d = 0;
totalt_s = 0;
for i = testnr
    test_name = tests(i).name;
    test_name = test_name(1:end-4);
    if( strcmp(statstr_d.(test_name),'OK') && strcmp(statstr_s.(test_name),'OK') )
       totalt_d = totalt_d + time_d.(test_name);
       totalt_s = totalt_s + time_s.(test_name);
       k = k+1;
    end
end
fprintf('Total time SEDUMI: %4.2f seconds, ECOS: %4.2f seconds (%4.2fx)\n', totalt_s, totalt_d, totalt_s/totalt_d);
