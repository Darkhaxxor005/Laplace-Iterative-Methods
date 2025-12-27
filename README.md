# 2D Laplace Equation – Iterative Methods

This project solves the 2D Laplace equation using iterative numerical methods and compares their convergence behavior.

---

## Methods Implemented
- Jacobi Method  
- Gauss–Seidel Method  
- Successive Over-Relaxation (SOR)

Finite difference discretization with Dirichlet boundary conditions is used.

---

## Project Structure

```text
Laplace-Iterative-Methods/
├── src/
│   ├── solvers.py
│   └── main.py
├── results/
│   └── plots/
├── README.md
├── requirements.txt 
```


---

## How to Run

Install dependencies:

```text
pip install -r requirements.txt
```

Run the program:

```text
python src/main.py
```

---

## User Inputs

During execution, the user provides:
- Grid size `N`
- Maximum iterations `max_iter`
- Convergence tolerance `tol`
- SOR relaxation factor `omega`

```text
The tol input structre means

1e-3 → 0.001 [ne-3 → 0.00n]
1e-6 → 0.000001 [ne-6 → 0.00000n]
1e-8 → 0.00000001 [ne-8 → 0.0000000n]
```

These parameters control accuracy and convergence behavior.

---

## Results

The program generates:
- Error vs iteration plot comparing all three methods  
  (`results/plots/error_comparison.png`)
- Heatmap of steady-state potential distribution  
  (`results/plots/potential_heatmap.png`)

---

## Observations

- Jacobi converges slowest  
- Gauss–Seidel converges faster  
- SOR converges fastest  

Results are consistent with theoretical expectations.

---

## Libraries Used

This project is implemented in Python and uses the following libraries:

### 1. NumPy
- Used for numerical computations
- Stores the 2D grid as arrays
- Used to compute error norms for convergence checking

### 2. Matplotlib
- Used for visualization
- Generates:
  - Error vs iteration convergence plots
  - Heatmap of the steady-state potential distribution

## Python Version
- Python 3.8 or higher is recommended
- I used Python 3.13.9

