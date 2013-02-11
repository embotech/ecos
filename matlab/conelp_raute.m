function z = conelp_raute(x,y,dims)
% Implements the "raute" operator for use in CONELP solver.
%
% z = conelp_raute(x,y,dims) returns z = x raute y
%
% NOTE: The solver and the text above are heavily based on the document
%
%  [1] L. Vandenberghe: "The CVXOPT linear and quadratic cone program 
%      solvers", March 20, 2010. 
%      [Online]: http://abel.ee.ucla.edu/cvxopt/documentation/coneprog.pdf
%  
%
% (c) Alexander Domahidi, IfA, ETH Zurich, 2012.

% LP cone
z = y(1:dims.l) ./ x(1:dims.l);

% SOC cone
for k = 1:length(dims.q)
    coneidx = dims.l+sum(dims.q(1:k-1))+1:dims.l+sum(dims.q(1:k));    
    xk = x(coneidx); yk = y(coneidx);
    x0 = xk(1); x1 = xk(2:end);    
    y0 = yk(1); y1 = yk(2:end);    
    xres = x0^2 - x1'*x1;
%     zk = 1/xkres*[xk0, -xk1'; -xk1, (xkres*eye(dims.q(k)-1) + xk1*xk1')/xk0]*yk;
    
    % fast linear code
    zeta = x1'*y1;
    z0 = x0*y0 - zeta;
    z1 = (zeta/x0 - y0).*x1 + (xres/x0).*y1;    
    z(coneidx,1) = [z0; z1] ./ xres;
end
