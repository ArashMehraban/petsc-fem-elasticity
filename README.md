# petsc-fem-elasticity
Petsc FEM code based on Brown Model

To run the code:
make all
mpiexec -n 2 ./main -polydegree 2 -dof 1 -f cube.exo -interpolate
