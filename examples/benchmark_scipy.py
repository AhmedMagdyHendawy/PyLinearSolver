from PyLinearSolver import IterativeSolvers
import numpy as np
from scipy.sparse import linalg
import time

# Positive Definite Matrix

# N= 2
n=2
A = np.random.rand(n,n) 
A = A * A.T 
b = np.random.rand(n)

start = time.time()
x = IterativeSolvers.cg(A,b)
print("Iterative Solvers CG Computational Time: %.7f"%(time.time()-start))
print("Iterative Solvers CG Solution: ",x)

start = time.time()
x = linalg.cg(A,b)
print("Scipy CG Computational Time: %.7f"%(time.time()-start))
print("Scipy CG Solution: ",x)

start = time.time()
x = IterativeSolvers.gmres(A+2,b)
print("Iterative Solvers CG Computational Time: %.7f"%(time.time()-start))
print("Iterative Solvers CG Solution: ",x)