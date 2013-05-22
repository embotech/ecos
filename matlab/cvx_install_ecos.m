makemex
[ ver, isoctave, fs, ps ] = cvx_version;

copyfile('cvx_ecos.m',strcat(fs,'/shims'));
mkdir(strcat(fs,'/ecos'));
copyfile(strcat('ecos.',ps),strcat(fs,'/ecos'));
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
