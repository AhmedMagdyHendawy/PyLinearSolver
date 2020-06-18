Installation
============

Dependencies
------------

To work with PyLinearSolver packages, you need to install the dependencies for them.

**Iterative Solvers**

    * Install `Julia <https://julialang.org/downloads/>`_
    * It is preferable to install a custom python version using `pyenv <https://github.com/pyenv/pyenv>`_ due to an issue `<https://pyjulia.readthedocs.io/en/latest/troubleshooting.html#ultimate-fix-build-your-own-python>`_ with PyJulia.


PyLinearSolver
--------------

Simply you need to call the following command::

    pip install PyLinearSolver

**NOTE** In case of working with Iterative Solvers, you need to install the dependencies for PyJulia through the following::

    import julia
    julia.install()
    