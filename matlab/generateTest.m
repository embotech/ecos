function generateTest(name,testdir,c,G,h,dims,A,b)
% Generates C test for ECOS from test data.
%
%     GenerateTest(name,testdir,c,G,h,dims,A,b) dumps the ECOS test problem
%     to the file <ecos>/test/testdir/name.h and creates a test such that
%     it can be included into the testsuite of ecostester.c
%
% (c) Alexander Domahidi, embotech GmbH, March 2014.


disp('Generating Test...');
txt = {'#include "ecos.h"'};
txt = [txt; '#include "minunit.h"'];
txt = [txt; sprintf('static char * test_%s(){',name)];
if( nargin > 6 )
    txt = [txt; cg_dump_conelpproblem(c,G,h,dims,A,b)];
elseif( nargin == 6 )
    txt = [txt; cg_dump_conelpproblem(c,G,h,dims)];
else
    error('generateTest requires at least 6 arguments. See help generateTest');
end
txt = [txt; ' '];
txt = [txt; 'pwork *mywork;'];
txt = [txt; 'idxint exitflag;'];
txt = [txt; ' '];
txt = [txt; '/* set up data */'];
txt = [txt; 'mywork = ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);'];
txt = [txt; 'if( mywork != NULL ){'];
txt = [txt; '/* solve */'];
txt = [txt; 'exitflag = ECOS_solve(mywork); }'];
txt = [txt; 'else exitflag = ECOS_FATAL;'];
txt = [txt; ' '];
txt = [txt; '/* clean up memory */'];
txt = [txt; 'ECOS_cleanup(mywork, 0);'];
txt = [txt; ' '];
txt = [txt; sprintf('mu_assert("%s: ECOS failed to produce outputflag OPTIMAL", exitflag == ECOS_OPTIMAL );',name)];
txt = [txt; 'return 0;'];
txt = [txt; '}'];

% save txt to name.h
fn = sprintf('%s/%s.h',testdir,name);
if( ~exist(testdir,'dir') )
    mkdir(testdir);
end
fprintf('Saving test to %s\n',fn);
cg_dumpfile(fn,txt);
