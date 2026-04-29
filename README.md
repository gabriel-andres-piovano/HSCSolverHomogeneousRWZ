# HSCSolverHomogeneousRWZ
Code for the computation of the homogeneous solutions of the Regge-Wheeler and Zerilli (RWZ) equations. The package includes an integrator that numerically solves the RWZ equation transformed in Hyperboloidal Slicing coordinates. All boundary conditions for the $`R^\text{in}`$ and $`R^\text{up}`$ solutions are implemented. The solver automatically switches to the asymptotic solutions near the horizon or $\infty$ for both $`R^\text{in}`$ and $`R^\text{up}`$.


Mathematica files
- ```HSCSolverHomogeneousRWZ.wl```: Mathematica package for the calculation of the numerical solutions to the Regge-Wheeler and Zerilli (RWZ) equations.
- ```Tutorial.nb```: short tutorial for the package

Installation
The file ```HSCSolverHomogeneousRWZ.wl``` must be placed at one of the paths shown in the variables $BaseDirectory (for system-wide installations) and $UserBaseDirectory (for single-user installations). You need to place it in the Applications/ subdirectory (or subfolder) of those returned by those variables. Alternatively, you can load the package using the command

```Get["/absolute_path_where_the_package_is_located/HSCSolverHomogeneousRWZ.wl"]; ```


Author:
- Gabriel Andres Piovano
