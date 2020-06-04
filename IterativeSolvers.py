from julia.api import Julia
from julia import Main

julia = Julia()

Main.using("IterativeSolvers")

# Wrapper function to Conjugate Gradients	
def cg(A,b,**kwargs):
    '''
    Arguments

        x: Initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.
    Keywords

        statevars::CGStateVariables: Has 3 arrays similar to x to hold intermediate results;
        initially_zero::Bool: If true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;
        Pl = Identity(): left preconditioner of the method. Should be symmetric, positive-definite like A;
        tol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol;
        maxiter::Int = size(A,2): maximum number of iterations;
        verbose::Bool = false: print method information;
        log::Bool = false: keep track of the residual norm in each iteration.
    Output

    if log is false

        x: approximated solution.
    if log is true

        x: approximated solution.
        ch: convergence history.
    ConvergenceHistory keys

        :tol => ::Real: stopping tolerance.
        :resnom => ::Vector: residual norm at each iteration.
    '''    
    return Main.cg(A,b,**kwargs)

# Wrapper function to Chebyshev iteration
def chebyshev(A,b,lamda_min,lamda_max,**kwargs):
    '''
    Arguments
        x: initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side;
        λmin::Real: lower bound for the real eigenvalues
        λmax::Real: upper bound for the real eigenvalues
    Keywords

        initially_zero::Bool = false: if true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;
        tol: tolerance for stopping condition |r_k| / |r_0| ≤ tol.
        maxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;
        Pl = Identity(): left preconditioner;
        log::Bool = false: keep track of the residual norm in each iteration;
        verbose::Bool = false: print convergence information during the iterations.
    Return values
        if log is false

            x: approximate solution.
        if log is true

            x: approximate solution;
            history: convergence history.
    '''
    return Main.chebyshev(A,b,lamda_min,lamda_max,**kwargs)

# Wrapper function to MINRES
def minres(A, b,**kwargs):
    '''
    Arguments
        x: initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.
    Keywords
        initially_zero::Bool = false: if true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;
        skew_hermitian::Bool = false: if true assumes that A is skew-symmetric or skew-Hermitian;
        tol: tolerance for stopping condition |r_k| / |r_0| ≤ tol. Note that the residual is computed only approximately;
        maxiter::Int = size(A, 2): maximum number of iterations;
        Pl: left preconditioner;
        Pr: right preconditioner;
        log::Bool = false: keep track of the residual norm in each iteration;
        verbose::Bool = false: print convergence information during the iterations.
    Return values
    if log is false
        x: approximate solution.
    if log is true
        x: approximate solution;
        history: convergence history.
    '''
    return Main.minres(A,b,**kwargs)

#Wrapper function to BiCGStab(l)
def bicgstabl(A, b, l, **kwargs):
    '''
    Arguments
        A: linear operator;
        b: right hand side (vector);
        l::Int = 2: Number of GMRES steps.

    Keywords
        max_mv_products::Int = size(A, 2): maximum number of matrix vector products.
        For BiCGStab(l) this is a less dubious term than "number of iterations";

        Pl = Identity(): left preconditioner of the method;
        tol::Real = sqrt(eps(real(eltype(b)))): tolerance for stopping condition |r_k| / |r_0| ≤ tol. Note that (1) the true residual norm is never computed during the iterations, only an approximation; and (2) if a preconditioner is given, the stopping condition is based on the preconditioned residual.
    Return values

    if log is false
        x: approximate solution.

    if log is true
        x: approximate solution;
        history: convergence history.
    '''
    return Main.bicgstabl(A,b,l,**kwargs)

#Wrapper function to IDR(s)
def idrs(A, b, s = 8, **kwargs):
    '''
    Arguments
        x: Initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.
        
    Keywords
        s::Integer = 8: dimension of the shadow space;
        tol: relative tolerance;
        maxiter::Int = size(A, 2): maximum number of iterations;
        log::Bool: keep track of the residual norm in each iteration;
        verbose::Bool: print convergence information during the iterations.
    Return values

    if log is false
        x: approximate solution.

    if log is true
    x: approximate solution;
        history: convergence history.
    '''
    return Main.idrs(A,b,s,**kwargs)

#Wrapper function to Restarted GMRES
def gmres(A, b, **kwargs):
    '''
    Arguments
        x: Initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.

    Keywords
        initially_zero::Bool: If true assumes that iszero(x) so that one matrix-vector product can be saved when computing the initial residual vector;
        tol: relative tolerance;
        restart::Int = min(20, size(A, 2)): restarts GMRES after specified number of iterations;
        maxiter::Int = size(A, 2): maximum number of inner iterations of GMRES;
        Pl: left preconditioner;
        Pr: right preconditioner;
        log::Bool: keep track of the residual norm in each iteration;
        verbose::Bool: print convergence information during the iterations.
    Return values

    if log is false
        x: approximate solution.
    if log is true
        x: approximate solution;
        history: convergence history.
    '''
    return Main.gmres(A,b,**kwargs)

# Wrapper function to LSMR - Least-squares minimal residual
def lsmr(A, b, **kwrags):
    '''
    Arguments

        x: Initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.
    Keywords
        λ::Number = 0: lambda.
        atol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).
        conlim::Number = 1e8: stopping tolerance. lsmr terminates if an estimate of cond(A) exceeds conlim. For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say). For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting
        atol = btol = conlim = zero, but the number of iterations may then be excessive.
        maxiter::Int = maximum(size(A)): maximum number of iterations.
        log::Bool: keep track of the residual norm in each iteration;
        verbose::Bool: print convergence information during the iterations.
    Return values

    if log is false

        x: approximated solution.
    if log is true

        x: approximated solution.
        ch: convergence history.
    ConvergenceHistory keys

        :atol => ::Real: atol stopping tolerance.
        :btol => ::Real: btol stopping tolerance.
        :ctol => ::Real: ctol stopping tolerance.
        :anorm => ::Real: anorm.
        :rnorm => ::Real: rnorm.
        :cnorm => ::Real: cnorm.
        :resnom => ::Vector: residual norm at each iteration.
    '''
    return Main.lsmr(A, b, **kwrags)

# Wrapper function to LSQR
def lsqr(A, b, **kwrags):
    '''
    Arguments

        x: Initial guess, will be updated in-place;
        A: linear operator;
        b: right-hand side.
    Keywords

        damp::Number = 0: damping parameter.
        atol::Number = 1e-6, btol::Number = 1e-6: stopping tolerances. If both are 1.0e-9 (say), the final residual norm should be accurate to about 9 digits. (The final x will usually have fewer correct digits, depending on cond(A) and the size of damp).
        conlim::Number = 1e8: stopping tolerance. lsmr terminates if an estimate of cond(A) exceeds conlim. For compatible systems Ax = b, conlim could be as large as 1.0e+12 (say). For least-squares problems, conlim should be less than 1.0e+8. Maximum precision can be obtained by setting atol = btol = conlim = zero, but the number of iterations may then be excessive.
        maxiter::Int = maximum(size(A)): maximum number of iterations.
        verbose::Bool = false: print method information.
        log::Bool = false: output an extra element of type ConvergenceHistory containing extra information of the method execution.
    Return values

    if log is false

        x: approximated solution.
    if log is true

        x: approximated solution.
        ch: convergence history.

    ConvergenceHistory keys
        :atol => ::Real: atol stopping tolerance.
        :btol => ::Real: btol stopping tolerance.
        :ctol => ::Real: ctol stopping tolerance.
        :anorm => ::Real: anorm.
        :rnorm => ::Real: rnorm.
        :cnorm => ::Real: cnorm.
        :resnom => ::Vector: residual norm at each iteration.
    '''
    return Main.lsqr(A, b, **kwrags)

## Stationary Methods

# Wrapper function to Jacobi
def jacobi(A, b, **kwargs):
    '''
        Keywords
            maxiter::Int = 10: maximum number of iterations.
    '''
    return Main.jacobi(A,b,**kwargs)

# Wrapper function to Gauss-Seidel
def gauss_seidel(A, b, **kwargs):
    '''
        Keywords
            maxiter::Int = 10: maximum number of iterations.
    '''
    return Main.gauss_seidel(A,b,**kwargs)

# Wrapper function to Successive over-relaxation (SOR)
def sor(A, b, w, **kwargs):
    '''
        Arguments
            w: relaxation parameter
        Keywords
            maxiter::Int = 10: maximum number of iterations.
    '''
    return Main.sor(A,b,w,**kwargs)

# Wrapper function to Symmetric successive over-relaxation (SSOR)
def ssor(A, b, w, **kwargs):
    '''
        Arguments
            w: relaxation parameter
        Keywords
            maxiter::Int = 10: maximum number of iterations.
    '''
    return Main.ssor(A,b,w,**kwargs)

