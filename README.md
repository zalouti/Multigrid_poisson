# Multigrid_poisson
> A Multigrid solver for Poisson equation implemented in the course Numeric Partial Differential Equation

## Files

- `gauss-seidel.m`: a linear equation solver using Gauss Seidel method.
- `grid_generator.m`: generate different levels of grids based on the finest spatial step and the number of grid level.
- `grid_interpolate.m`: interpolate grid points for fine mesh based on last level of coarse mesh.
- `matrix_assemble.m`: discrete Poisson equation in the area of $[0, 1]\times [0, 1]$
- `relax.m`: relax the equation on the current level of mesh
- `restrict.m`: transfer residual to the coarse mesh
- `main.m`: define parameters and equation, implement multigrid method and plot the final result

## Figure

![V-cycle convergence speed](<convergence speed.png>)
