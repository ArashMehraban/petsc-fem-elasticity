#ifndef __FE__H
#define __FE__H

#include<petscdt.h>
/* includes the header for PetscDTGaussQuadrature
 that calculates the quadrature points and weights*/

 typedef enum {Q1, Q2} RestricMode;

 typedef struct FE_private *FE;

 struct FE_private
{
  MPI_Comm comm;
  PetscInt polydegree; //Finite Element polynomial degree
  PetscInt dof;        //degrees of freedom at each node
  PetscInt      *connQ1,
                *connQ2,
                sz_connQ1,
                sz_connQ2,
                sz_perm_idx_Q1,
                sz_perm_idx_Q2;
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
