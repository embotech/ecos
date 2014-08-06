#!/usr/bin/env python
from qcml import QCML


if __name__ == '__main__':
    test_problems = [
    ("""variable x(3)
        minimize norm(x)
        x >= 0""", "norm"),
    ("""variable x(3)
        minimize quad_over_lin(x,1)
        x >= 0""", "quad_over_lin"),
    ("""variable x(3)
        minimize square(norm(x))
        x >= 0""", "sq_norm"),
    ("""variable x(3)
        minimize sum(square(x))
        x >= 0""", "sum_sq"),
    ("""variable x
        minimize inv_pos(x)""", "inv_pos")
    ]

    p = QCML(debug=True)
    def solve(s):
        print s[1]
        p.parse(s[0])
        p.canonicalize()
        p.codegen("C")
        p.save(s[1])
        # socp_data = p.prob2socp(params=locals())
        # import ecos
        # sol = ecos.solve(**socp_data)
        # return sol

    print map(solve, test_problems)
