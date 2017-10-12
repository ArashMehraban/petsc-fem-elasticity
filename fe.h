#ifndef __FE__H
#define __FE__H

#include<petscdt.h>
/* includes the header for PetscDTGaussQuadrature
 that calculates the quadrature points and weights*/

 typedef struct FE_private *FE;

 struct FE_private
{
  MPI_Comm comm;
  PetscInt polydegree;   //Finite Element polynomial degree
  PetscInt addquadpts;   // Number of additonal quadrature points
  PetscInt dof;          //degrees of freedom at each node
  struct{
    PetscReal *B;
    PetscReal *D;
    PetscReal *D_tilda;
    PetscReal *x;
    PetscReal *w;
    PetscReal *w3;
  }ref;
};

PetscErrorCode FEbasisEval(FE fe, PetscReal q, PetscReal B[], PetscReal D[]);
PetscErrorCode FESetup(FE fe);
PetscErrorCode dmFEcreate(PetscInt dof, PetscInt polydegree, PetscInt addquadpts, FE *fe);

#endif //end of __FE__H
