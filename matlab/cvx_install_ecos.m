makemex
fs = cvx_where;

copyfile('cvx_solver_shim.m',strcat(fs,filesep, 'shims', filesep, 'cvx_ecos.m'));
mkdir(strcat(fs,filesep,'ecos'));
copyfile(strcat('ecos.',mexext),strcat(fs,filesep,'ecos'));
cvx_setup

% mini test script
%{
cvx_begin
cvx_solver 'ecos'
variable x
minimize(x)
x>=1
cvx_end
%}
