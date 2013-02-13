function L = conelp_lowrankL(p,q)
% Returns the matrix
%
%                    [ 1                                 ]
%                    [ p2*q1   1                         ]
%                    [ p3*q1  p3*q2    1                 ]
%          L(p,q) =  [   .      .           .            ]
%                    [   .      .               .        ]
%                    [   .      .                   .    ]
%                    [ pn*q1  pn*q2   . . . pn*qn-1   1  ]
%
% This function is for testing purposes only. Note that L(p,q) should never
% be explicitly formed during the factorization/solve process, since this
% is the whole point of doing low-rank updates.
%
% (c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2013.
% Email: domahidi@control.ee.ethz.ch


n = length(p);

L = eye(n);
for j = 1:n-1
    L(j+1:end,j) = q(j).*p(j+1:end);
end
