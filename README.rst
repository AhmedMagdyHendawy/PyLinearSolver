.. image:: https://readthedocs.org/projects/pylinearsolver/badge/?version=latest
    :target: https://pylinearsolver.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

About the Package
=================

PyLinearSolver is a python iterface for many well-known linear solvers packages written in different languages like C, C++, julia, etc...
The package uses the existing implementation from the original source and build a wrapper layer for a python development. 

**NOTE** All rights are preserved to the authers and developers of the original packages.


Installation
============

Dependencies
------------

To work with PyLinearSolver packages, you need to install the dependencies for them.

**Iterative Solvers**

    * Install `Julia <https://julialang.org/downloads/>`_


PyLinearSolver
--------------

Simply you need to call this command::

    pip install PyLinearSolver

**NOTE** In case of working with Iterative Solvers, you need to install the dependencies for PyJulia through the following::

    import julia
    julia.install()
    
    
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