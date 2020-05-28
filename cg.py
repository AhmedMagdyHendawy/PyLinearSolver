# from julia import Main
from julia.api import Julia
import numpy as np 
julia = Julia(compiled_modules=False)
from julia import Main
Main.using("IterativeSolvers")


def cg(A,b,verbose=False,log=False):    
    return Main.cg(A,b,verbose=verbose,log=log)

n=10
A = np.random.rand(n,n)
b= np.random.rand(n)

A = A + A.T +2*n*np.identity(n)

x, ch=cg(A,b,verbose=True, log=True)

