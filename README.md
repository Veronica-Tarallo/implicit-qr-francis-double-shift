# Implicitly Shifted QR Algorithm (Francis Double Shift)

A from-scratch MATLAB implementation of the **Implicitly Shifted QR Algorithm with Francis Double Shift** to compute the full eigenvalue spectrum of generic real non-symmetric matrices in $\mathcal{O}(n^3)$ operations.

## Overview & Mathematical Formulation

Computing eigenvalues of non-symmetric matrices requires handling complex-conjugate pairs efficiently while avoiding numerical instability. This solver implements the standard pipeline used in modern numerical linear algebra libraries:

1. **Hessenberg Reduction:** Reduces a dense matrix $A \in \mathbb{R}^{n \times n}$ to upper Hessenberg form $H$ via orthogonal similarity transformations (Householder reflections), lowering per-iteration cost from $\mathcal{O}(n^3)$ to $\mathcal{O}(n^2)$.
2. **Francis Double-Shift Step (Implicit Bulge Chasing):** Uses the trace and determinant of the trailing $2 \times 2$ submatrix to implicitly apply two QR steps simultaneously. By introducing a small perturbation and chasing the resulting $3 \times 1$ bulge down the subdiagonal with Householder reflectors, complex-conjugate eigenvalue pairs are computed entirely in **real arithmetic**.
3. **Deflation Criteria:** Monitors subdiagonal entries against a relative tolerance of $10^{-6}$ to decouple converged $1 \times 1$ real eigenvalues and $2 \times 2$ real Schur blocks representing complex-conjugate pairs.

## Numerical Validation & Benchmarks

The implementation was validated against MATLAB's built-in `eig` function across randomized test matrices:

* **Accuracy:** Achieved a maximum absolute eigenvalue error of $\sim 5.4 \times 10^{-8}$.
* **Convergence Speed:** Empirically confirmed $\mathcal{O}(1)$ average iterations per deflation (~2.64 iterations per eigenvalue/pair).
* **Asymptotic Complexity:** Benchmarked runtime scaling across matrix sizes $n = 10$ to $n = 90$, verifying the expected $\mathcal{O}(n^3)$ computational complexity.
* **Stagnation Stress-Test ($n = 400$):** Tested on larger non-symmetric matrices up to $n = 400$ to isolate numerical stagnation phenomena, documenting how production libraries (such as LAPACK) inject *exceptional shifts* to break symmetry loops and guarantee global convergence.

## Repository Structure

* `main.m` - Main driver script running convergence plots, accuracy checks, asymptotic complexity benchmarks, and stagnation stress-tests.
* `implicitQRAlg.m` - Core solver implementing upper Hessenberg reduction, Francis Double Shift, bulge chasing, and $1 \times 1$ / $2 \times 2$ deflation.
* `householder_vector.m` - Helper function computing normalized Householder reflectors.

## How to Run

1. Clone or download this repository.
2. Open MATLAB and navigate to the project folder.
3. Run `main` in the MATLAB Command Window to execute the numerical experiments and display benchmark results.
