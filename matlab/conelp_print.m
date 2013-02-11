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
        case 'ldlsparse', disp('Using LDLSPARSE as linear solver on unreduced system');        
        case 'ldlsparse2', disp('Using LDLSPARSE as linear solver on unreduced system, avoiding W^2');     
        case 'backslash', disp('Using MATLAB''s backslash as linear system solver on unreduced system');
        case 'backslash2', disp('Using MATLAB''s backslash as linear system solver on reduced system');
        case 'backslash3', disp('Using MATLAB''s backslash as linear system solver on full system, avoiding W^2');    
    end
    fprintf('\n');
    if( ismac )
        fprintf('It     pcost       dcost     gap    pres   dres    k/t    mu    stepA  step\n');    
    else
        fprintf('It     pcost         dcost      gap     pres    dres    k/t      mu    stepA  step\n');    
    end
    fprintf('%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e   N/A    N/A \n',info.iter,info.pcost,info.dcost,info.gap,info.pres,info.dres,info.kapovert,info.mu);
else    
    fprintf('%2d  %+5.3e  %+5.3e  %+2.0e  %2.0e  %2.0e  %2.0e  %2.0e  %5.3f  %5.4f\n',info.iter,info.pcost,info.dcost,info.gap,info.pres,info.dres,info.kapovert,info.mu,info.step_aff,info.step);
end

