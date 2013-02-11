function alpha = conelp_linesearch_soc(alpha, s, ds)
% Returns the "best" intersection of step sizes for a second order cone.

a = ds(1)^2 - ds(2:end)'*ds(2:end);
b = 2*(ds(1)*s(1) - ds(2:end)'*s(2:end));
c = s(1)^2 - s(2:end)'*s(2:end);
d = b^2 - 4*a*c;

assert(a ~= 0, 'Divisor is zero in line search on s');
if( a > 0 && c > 0 )
     newa = [-inf, +inf];
elseif( a < 0 && c < 0 )
        error('Infeasible search direction - no progress possible.');
else
    newa(1,1) = (-b + sqrt(d))/(2*a);
    newa(1,2) = (-b - sqrt(d))/(2*a);
end
assert(all(isreal(newa)),'Complex line search parameter on s');
%newa( newa < 0 ) = 0;
%assert(sum(newa > 0) >= 1,'No progress possible - maximum step size is 0') );
alpha = min([alpha, newa( newa > 0 )]);
if(alpha < 0 )
    error('Line search exited with negative step length');
end