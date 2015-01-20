
%Load the lpnetlib problems and test the lp solver
	
  %The first 10 lps in standard form from LP netlib
  standard_form_problems = [ 594   596   597   598   599   600   601   602   603   604  ];
  %Structure for the results
  results = {{'Problem name','Iter','Flag','Norm Residual','Complementarity','Conic Infeasiblity'}}

  %for problem_index = 1:length(standard_form_problems)

  for problem_index = 1:10
    
    %Get the problem form UFGET
    problem_uf_ix = standard_form_problems(problem_index);
    P = UFget(problem_uf_ix);
      
	fprintf('Problem name %s\n',P.name)
    A = P.A; 
    b = P.b;
    c = P.aux.c;
	
    [p,n] = size(A);
	m = n;
	G     = -speye(n);
	h     = zeros(n,1);
	%Set the dimensions
	dims.l = n;
	[x,y,info,s,z] = ecos(c,G,h,dims,A,b);
	%Check the residuals
	pres = norm([A*x-b;G*x+s-h]);
	dres = norm(A'*y+G'*z+c);
	lin_res = sqrt(pres^2+dres^2)
	comp = s'*z;
	conic_infeas = min(min(s),min(z));

	fprintf('Residuals: primal %3.3e, dual %3.3e complementarity %3.3e conic infeas %3.3e \n',pres,dres,comp,conic_infeas)

	results = {results{:},{P.name,info.iter,info.exitflag,lin_res,comp,conic_infeas}}; 
end

res = results{1};
fprintf('%20s %s %s %s %s %s\n',res{1},res{2},res{3},res{4},res{5},res{6});
for(j=2:length(results))
    res = results{j};
    fprintf('%20s, %3i, %3i, %2.1e, %2.1e %2.1e\n',res{1},res{2},res{3},res{4},res{5},res{6});
end
