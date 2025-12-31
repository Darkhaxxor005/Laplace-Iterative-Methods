# src/solvers.py
""""    
    Parameters
    ----------
    u : Initial grid with boundary conditions applied
    max_iter : Maximum number of iterations
    tol : Convergence tolerance
    omega : Relaxation factor (1 < omega < 2)
    
    Returns
    -------
    u : Converged solution
    errors : Error at each iteration
    """
import numpy as np


def apply_boundary_conditions(u):
    """
    Apply Dirichlet boundary conditions to the grid.

    Top boundary    : 100
    Bottom boundary : 0
    Left boundary   : 0
    Right boundary  : 0
    """
    u[0, :] = 100  # Defines boundary values.
    u[-1, :] = 0
    u[:, 0] = 0
    u[:, -1] = 0
    return u


def jacobi(u, max_iter, tol):
    
    # Solve 2D Laplace equation using the Jacobi iterative method.

    u_new = u.copy()
    errors = []

    for _ in range(max_iter):
        for i in range(1, u.shape[0] - 1):
            for j in range(1, u.shape[1] - 1):
                u_new[i, j] = 0.25 * (
                    u[i + 1, j] + u[i - 1, j] +      # Finite Difference Discretization Equation
                    u[i, j + 1] + u[i, j - 1]
                )

        error = np.linalg.norm(u_new - u)   # The error is defined as the L2 norm of the difference between successive solution iterates.
        errors.append(error)

        if error < tol:
            break  # Breaks the loop when the error is smaller than the error threshold.

        u = apply_boundary_conditions(u_new.copy())

    return u, errors   # Returns u and errors to the main function.


def gauss_seidel(u, max_iter, tol):

    # Solve 2D Laplace equation using the Gauss-Seidel method.

    errors = []

    for _ in range(max_iter):
        u_old = u.copy()

        for i in range(1, u.shape[0] - 1):
            for j in range(1, u.shape[1] - 1):
                u[i, j] = 0.25 * (
                    u[i + 1, j] + u[i - 1, j] +     # Finite Difference Discretization Equation
                    u[i, j + 1] + u[i, j - 1]
                )

        error = np.linalg.norm(u - u_old)    # The error is defined as the L2 norm of the difference between successive solution iterates.
        errors.append(error)

        if error < tol:
            break    # Breaks the loop when the error is smaller than the error threshold.

        u = apply_boundary_conditions(u)

    return u, errors  # Returns u and errors to the main function.


def sor(u, omega, max_iter, tol):

    # Solve 2D Laplace equation using the SOR method.

    errors = []

    for _ in range(max_iter):
        u_old = u.copy()

        for i in range(1, u.shape[0] - 1):
            for j in range(1, u.shape[1] - 1):
                u[i, j] = (1 - omega) * u[i, j] + omega * 0.25 * (
                    u[i + 1, j] + u[i - 1, j] +       # Finite Difference Discretization Equation
                    u[i, j + 1] + u[i, j - 1]
                )

        error = np.linalg.norm(u - u_old)     # The error is defined as the L2 norm of the difference between successive solution iterates.
        errors.append(error)

        if error < tol:
            break   # Breaks the loop when the error is smaller than the error threshold.

        u = apply_boundary_conditions(u)

    return u, errors  # Returns u and errors to the main function.
