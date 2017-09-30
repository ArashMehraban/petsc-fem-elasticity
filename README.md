# petsc-fem-elasticity
Petsc FEM code based on Brown Model

To run the code:
make all
mpiexec -n 2 ./main -f cube8.exo -ne 4 -dof 3 -polydegree 2 
