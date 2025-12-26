# src/main.py
import numpy as np
import matplotlib.pyplot as plt
from solvers import jacobi, gauss_seidel, sor, apply_boundary_conditions

# --------------------------------
# User Input
# --------------------------------
print("2D Laplace Equation Solver using Iterative Methods")

N = int(input("Enter grid size N (e.g. 50): "))
max_iter = int(input("Enter maximum iterations (e.g. 5000): "))
tol = float(input("Enter tolerance (e.g. 1e-6): "))
omega = float(input("Enter SOR relaxation factor omega (e.g. 1.7): "))

# --------------------------------
# Initial Grid
# --------------------------------
u0 = np.zeros((N, N))
u0 = apply_boundary_conditions(u0)

# --------------------------------
# Solve using Jacobi, Gauss–Seidel, and SOR.
# --------------------------------
u_jacobi, err_j = jacobi(u0.copy(), max_iter, tol)
u_gs, err_gs = gauss_seidel(u0.copy(), max_iter, tol)
u_sor, err_sor = sor(u0.copy(), omega, max_iter, tol)

# --------------------------------
# Plotting Results
# --------------------------------
plt.figure()
plt.plot(err_j, label="Jacobi")
plt.plot(err_gs, label="Gauss-Seidel")
plt.plot(err_sor, label="SOR")
plt.yscale("log")
plt.xlabel("Iteration")
plt.ylabel("Error")
plt.legend()
plt.title("Convergence Comparison")
plt.savefig("results/plots/error_comparison.png")

plt.figure()
plt.imshow(u_sor, cmap="hot")
plt.colorbar()
plt.title("Steady State Potential (SOR)")
plt.savefig("results/plots/potential_heatmap.png")

plt.show()
