% New CVX tests ought to be *scripts* which name the variables to be
% compared to "ecos_NAME" and "true_NAME", where "true_NAME" is computed
% from SeDuMi or known ahead of time.
%
% This script will automatically compare "ecos_NAME" to "true_NAME".
%
% Optionally, if 'time_ecos' and 'time_true' are defined, the script will
% display this information after a notification of success.
%
clear all;
tests = dir('test_*.m');

TOL = 1e-4;

disp(' ')
for i = 1:length(tests),
    cvx_clear
    clear ecos_* true_* time_*;
    % remove .m from end of script name
    script = tests(i).name(1:end-2);
    disp(['TESTING ''' script '''']);
    evalc( script );
    to_compare = who('ecos_*');
    timings = who('time_*');
    if ~isempty(timings)
        ecos_time = time_ecos;
        true_time = time_true;
        timings = sprintf('ecos (%0.4f), default (%0.4f)', ecos_time, true_time);
        disp(['  TIMING: ' timings]);
    end
    for j = 1:length(to_compare),
        ecos_str = to_compare{j};
        true_str = strrep(ecos_str, 'ecos', 'true');
        ecos_val = eval(ecos_str);
        true_val = eval(true_str);
        n = length(ecos_val);
        err = norm(ecos_val - true_val);
        try
            assert(err <= TOL * sqrt(n))
            disp('  SUCCESS')
        catch e
            err_msg = sprintf('  TEST FAIL: ''%s'' and ''%s'' do not agree in ''%s''', ecos_str, true_str, script);
            disp(err_msg)
            disp(['    norm error of      : ' num2str(err)]);
            disp(['    expected less than : ' num2str(TOL * sqrt(n))]);
        end
    end
        
end