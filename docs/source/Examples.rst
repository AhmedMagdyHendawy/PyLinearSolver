Examples
========

Iterative Solvers
-----------------

An example for Iterative Solver Wrapper package in PyLinearSolver::

    from PyLinearSolver import IterativeSolvers
    import numpy as np

    n=10
    A = np.random.rand(n,n)
    b= np.random.rand(n)
    x= np.zeros(n)
    A = A + A.T +2*n*np.identity(n)

    x, ch=IterativeSolvers.cg(A,b,verbose=True, log=True)