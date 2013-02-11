function makemex(what)

if( nargin == 0 )
    what = {'all'};
elseif( isempty(strfind(what, 'all')) && ...
        isempty(strfind(what, 'ldl')) && ...
        isempty(strfind(what, 'amd')) && ...
        isempty(strfind(what, 'clean')) && ...
        isempty(strfind(what, 'ecos')) && ...
        isempty(strfind(what, 'ecosmex')) )
    fprintf('No rule to make target "%s", exiting.\n', what);
end


if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims -DDLONG -DLDL_LONG';    
else
    d = '-DDLONG -DLDL_LONG';    
end


% ldl solver
if( any(strcmpi(what,'ldl')) || any(strcmpi(what,'all')) )
    fprintf('Compiling LDL...');
    cmd = sprintf('mex -c -O %s -I../external/ldl/include -I../external/SuiteSparse_config ../external/ldl/src/ldl.c', d);
    eval(cmd);
    fprintf('\t\t\t[done]\n');
end


% amd (approximate minimum degree ordering)
if( any(strcmpi(what,'amd')) || any(strcmpi(what,'all')) )
    fprintf('Compiling AMD...');
    i = sprintf ('-I../external/amd/include -I../external/SuiteSparse_config') ;
    cmd = sprintf ('mex -c -O %s -DMATLAB_MEX_FILE %s', d, i) ;
    files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
             'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
             'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
    for i = 1 : length (files)
        cmd = sprintf ('%s ../external/amd/src/%s.c', cmd, files {i}) ;
    end
    eval(cmd);
    fprintf('\t\t\t[done]\n');
end


% ecos (solver)
if( any(strcmpi(what,'ecos')) || any(strcmpi(what,'all')) ) 
    fprintf('Compiling ecos...');
    i = sprintf('-I../include -I../external/SuiteSparse_config -I../external/ldl/include -I../external/amd/include');
    cmd = sprintf('mex -c -O %s -DMATLAB_MEX_FILE %s', d, i);
    files = {'ecos', 'kkt', 'cone', 'spla', 'timer', 'preproc', 'splamm'};
    for i = 1 : length (files)
        cmd = sprintf ('%s ../src/%s.c', cmd, files {i}) ;
    end
    eval(cmd);    
    fprintf('\t\t\t[done]\n');
end


% ecos_mex
if( any(strcmpi(what,'ecosmex')) || any(strcmpi(what,'all')) )
    fprintf('Compiling ecos_mex...');
    cmd = sprintf('mex -c -O -DMATLAB_MEX_FILE -DMEXARGMUENTCHECKS %s -I../include -I../external/SuiteSparse_config -I../external/ldl/include -I../external/amd/include ecos_mex.c', d);
    eval(cmd);
    fprintf('\t\t\t[done]\n');
    fprintf('Linking...     ');
    clear ecos
    if( ispc )
        cmd = sprintf('mex %s amd_1.obj amd_2.obj amd_aat.obj amd_control.obj amd_defaults.obj amd_dump.obj amd_global.obj amd_info.obj amd_order.obj amd_post_tree.obj amd_postorder.obj amd_preprocess.obj amd_valid.obj ldl.obj kkt.obj preproc.obj spla.obj cone.obj ecos.obj timer.obj splamm.obj ecos_mex.obj -output "ecos"', d);
        eval(cmd);    
    elseif( ismac )
        cmd = sprintf('mex %s -lm amd_1.o   amd_2.o   amd_aat.o   amd_control.o   amd_defaults.o   amd_dump.o   amd_global.o   amd_info.o   amd_order.o   amd_post_tree.o   amd_postorder.o   amd_preprocess.o   amd_valid.o     ldl.o   kkt.o   preproc.o   spla.o   cone.o   ecos.o timer.o   splamm.o   ecos_mex.o   -output "ecos"', d);
        eval(cmd);
    elseif( isunix )
        cmd = sprintf('mex %s -lm -lrt amd_1.o   amd_2.o   amd_aat.o   amd_control.o   amd_defaults.o   amd_dump.o   amd_global.o   amd_info.o   amd_order.o   amd_post_tree.o   amd_postorder.o   amd_preprocess.o   amd_valid.o     ldl.o   kkt.o   preproc.o   spla.o   cone.o   ecos.o timer.o   splamm.o   ecos_mex.o   -output "ecos"', d);
        eval(cmd);
    end
    fprintf('\t\t\t\t[done]\n');
        
%     fprintf('Copying MEX file...');
%     clear ecos
%     copyfile(['ecos.',mexext], ['../ecos.',mexext], 'f');
%     copyfile( 'ecos.m', '../ecos.m','f');
%     fprintf('\t\t\t[done]\n');
    disp('ecos successfully compiled. Happy solving!');
end

  
% clean
if( any(strcmpi(what,'clean')) || any(strcmpi(what,'purge')) )
    fprintf('Cleaning up object files...  ');
    if( ispc ), delete('*.obj'); end
    if( isunix), delete('*.o'); end
    fprintf('\t[done]\n');
end

% purge
if( any(strcmpi(what,'purge')) )
    fprintf('Deleting mex file...  ');
    if( ispc ), delete(['ecos.',mexext]); end
    if( isunix), delete(['ecos.',mexext]); end
    fprintf('\t\t[done]\n');
end
