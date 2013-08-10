#/usr/bin/python
from scoop import QCML
import argparse

parser = argparse.ArgumentParser(description="generate code for testing")
parser.add_argument('-m', nargs=1, default=[1],type=int,dest='m', help="problem dimension")
parser.add_argument('-n', nargs=1, default=[1],type=int,dest='n', help="size of variable")
parser.add_argument('-q', nargs=1, default=[None],type=int,dest='q', help="cone size")
parser.add_argument('--cvx', help="solve the problem using CVX", action="store_true")

group = parser.add_mutually_exclusive_group()
group.add_argument("--svm", help="solve the SVM problem", action="store_true")
group.add_argument("--portfolio", help="solve the portfolio problem", action="store_true")

args = parser.parse_args()

p = QCML()
m = args.m[0]
n = args.n[0]
q = args.q[0]

if args.svm:
    
    p.parse("""
        dimensions m n
        variable a(n)
        variable b
        parameter X(m,n)      # positive samples
        parameter Y(m,n)      # negative samples
        parameter gamma positive
    
        minimize (norm(a) + gamma*sum(pos(1 - X*a + b) + pos(1 + Y*a - b)))
    """)
    
    p.rewrite()
    if args.cvx:
        p.codegen("cvx")
    else:
        p.codegen("matlab",cone_size=q,m=m,n=n)
    p.prettyprint()

elif args.portfolio:
    p.parse("""
        dimensions m n
        
        variable x(n)
        parameter mu(n)
        parameter gamma positive
        parameter F(n,m)
        parameter D(n,n)
        maximize (mu'*x - gamma*(square(norm(F'*x)) + square(norm(D*x))))
            sum(x) == 1
            x >= 0
    """)
    
    p.rewrite()
    if args.cvx:
        p.codegen("cvx")
    else:
        p.codegen("matlab",cone_size=q,m=m,n=n)
    p.prettyprint()
    

else:
    print "you must decide whether you want to solve the SVM or Portfolio problem."
