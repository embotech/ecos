function conelp_print(info)
% Printing utility function or CONELP solver.
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

if( info.iter == 0 )
    % first iteration, print header as well
    fprintf('\n    ***************************************************************************\n');
    fprintf('    * This is CONELP, a MATLAB interior point solver for convex cone programs *\n');
    fprintf('    *                                                                         *\n');
    fprintf('    * NOTE: The solver is heavily based on L. Vandenberghe''s                  *\n');
    fprintf('    * "The CVXOPT linear and quadratic cone program solvers", March 20, 2010. *\n');
    fprintf('    * [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf     *\n');
    fprintf('    *                                                                         *\n');
    fprintf('    * (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2012. *\n');
    fprintf('    *     Email: domahidi@control.ee.ethz.ch                                  *\n');
    fprintf('    ***************************************************************************\n');
    fprintf('\n');
    switch( lower(info.LINSOLVER) )
        case 'ldlsparse',    disp('Search direction computation: LDLSPARSE with sparse SOC scalings');  
        case 'rank1updates', disp('Search direction computation: LDLSPARSE + rank 1 updates');  
        case 'backslash',    disp('Search direction computation: MATLAB''s LDL + backslash');
    end
    fprintf('\n');
    if( ismac )
        fprintf('It     pcost       dcost     gap    pres   dres    k/t    mu    stepA  step     IR \n');    
    else
        fprintf('It     pcost         dcost      gap     pres    dres    k/t      mu    stepA  step    IR\n');    
    end
    fprintf('%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e   N/A    N/A     N/A\n',info.iter,info.pcost,info.dcost,info.gap,info.pres,info.dres,info.kapovert,info.mu);
else    
    fprintf('%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e  %5.3f  %5.4f  %d %d %d\n',info.iter,info.pcost,info.dcost,info.gap,info.pres,info.dres,info.kapovert,info.mu,info.step_aff,info.step,info.nitref1,info.nitref2,info.nitref3);
end

