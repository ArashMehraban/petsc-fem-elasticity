#ifndef __USER__H
#define __USER__H

#include<petscdmplex.h>
/*for dm and dmplex objects*/
#include "fe.h"
#include "math.h"
#include "topology.h"


typedef struct{
  char          filename[PETSC_MAX_PATH_LEN]; //optional exodusII file
  PetscBool	    interpolate;                  //Generate intermediate mesh points: nodes on faces, edges in addition to vertices and element
  PetscInt      polydegree;                   //polynomial degree (1 less than number of interpolating points in 1D)
  PetscInt      addquadpts;                   // Number of additonal quadrature points
  PetscInt      dof;                          //number of degrees of freedom at each node, e.g. For 3D linear elasticity: at each node we have u1, u2 and u3
  PetscInt      ne;                           //number of elements to extract per iteration
}AppCtx;

PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx *userOptions);
PetscErrorCode dmCreate(MPI_Comm comm, AppCtx user, DM *dm);
/* Use drawOneElem after calling dmCreate.
   It is best to use this function with one element to see its connectivity */
PetscErrorCode drawOneElem(DM dm,AppCtx user);

#endif //end of __USER__H
