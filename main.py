import IterativeSolvers
import numpy as np

n=10
A = np.random.rand(n,n)
b= np.random.rand(n)
x= np.zeros(n)
A = A + A.T +2*n*np.identity(n)

x, ch=IterativeSolvers.cg(A,b,verbose=True, log=True)
print(x)

x, ch=IterativeSolvers.chebyshev(A,b,0.1,0.2,verbose=True, log=True)
print(x)