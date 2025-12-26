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

pip install -r requirements.txt

Run the program:

python src/main.py

---

## User Inputs

During execution, the user provides:
- Grid size `N`
- Maximum iterations `max_iter`
- Convergence tolerance `tol`
- SOR relaxation factor `omega`

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
