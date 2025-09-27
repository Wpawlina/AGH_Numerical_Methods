AGH UST Course of Date Structures and Algorithms
=============================================== 
This repository contains the materials for the course of Numerical methods at the AGH University of Science and Technology,
algorithms learned in this course and exercises solved during the course.


## Algorithms in this repository

- Lab 1 — Numerical differentiation (finite differences)
	- Forward-difference derivative f'(x) ≈ (f(x+h) − f(x)) / h, with precision study for float/double/long double
	- File: `lab1/pawlina_solution.cpp`

- Lab 2 — Polynomial interpolation
	- Lagrange interpolation (uniform and Chebyshev nodes)
	- Newton interpolation via divided differences (uniform and Chebyshev nodes)
	- Error estimation and plotting utilities
	- File: `lab2/solution.py`

- Lab 3 — Hermite interpolation
	- Hermite polynomial interpolation using repeated nodes and divided differences (uniform and Chebyshev nodes)
	- Error estimation and plotting
	- File: `lab3/solution.py`

- Lab 4 — Spline interpolation and tridiagonal solver
	- Natural/clamped cubic splines (build and solve tridiagonal system)
	- Quadratic splines
	- Thomas algorithm (tridiagonal solver) used in spline setup
	- File: `lab4/soulution.py`

- Lab 5 — Least-squares approximation and linear solver
	- Polynomial least-squares approximation (normal equations)
	- Gaussian elimination with partial pivoting (custom implementation)
	- Error estimation and plotting
	- File: `lab5/solution.py`

- Lab 6 — Fourier approximation
	- Discrete Fourier series approximation on an interval (compute a0, an, bn from sampled data)
	- Error estimation and plotting
	- File: `lab6/solution.py`

- Lab 7 — Nonlinear equation solvers
	- Newton–Raphson method (with residual- or step-based stopping)
	- Secant method (with residual- or step-based stopping)
	- File: `lab7/solution.py`

- Lab 8 — Linear systems (C++ and Python)
	- Thomas algorithm for tridiagonal systems (C++: float and double precision)
	- Gaussian elimination with partial pivoting (C++: float and double precision)
	- Jacobi iterative method (Python)
	- Matrix generators for test systems; timing and accuracy checks
	- Files:
		- `lab8/solution_float.cpp`
		- `lab8/solution_double.cpp`
		- `lab8/solution.py`

- Lab 9 — Linear systems analysis and solvers (Python)
	- Matrix generators (3 patterns), condition number analysis, spectral radius helper
	- Gaussian elimination (dense), Thomas algorithm (tridiagonal), Jacobi method (iterative)
	- File: `lab9/solution.py`

- Lab 10 — ODE initial value problems
	- Explicit Euler method
	- Runge–Kutta 4th-order (RK4) method
	- Error comparison and data export
	- File: `lab10/solution.py`


