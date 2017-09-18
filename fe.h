#ifndef __FE__H
#define __FE__H

#include<petscdt.h>
/* includes the header for PetscDTGaussQuadrature
 that calculates the quadrature points and weights*/

 typedef struct FE_private *FE;

 struct FE_private
{
  MPI_Comm comm;
  PetscInt polydegree; //Finite Element polynomial degree
  PetscInt dof;        //number of degrees of freedom per vertex
  struct{
    PetscReal *B;
    PetscReal *D;
    PetscReal *x;
    PetscReal *w;
    PetscReal *w3;
  }ref;
};

PetscErrorCode FEbasisEval(FE fe, PetscReal q, PetscReal B[], PetscReal D[]);
PetscErrorCode FESetup(FE fe);


#endif //end of __FE__H
