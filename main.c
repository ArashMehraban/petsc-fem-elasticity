static char help[] = "fem-elasticity";

#include<petscsnes.h>
/*petscsnes.h also includes the following headers:
  1) petscsys.h
  2) petscmat.h
  3) petscis.h
  4) petscviewer.h
  5) petscksp.h
  6) petscpc.h
*/
#include "fe-align.h"
#include "pointwise.h"
#include "fe.h"
#include "user.h"
#include "exact.h"

//PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
//PetscErrorCode FormJacobian(SNES,Vec, Mat,Mat, void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  //PetscMPIInt       rank;
  PetscErrorCode    ierr;
  DM                dm;
  AppCtx	          user; /*user-defined work context*/
  Vec               exactSol;
  FE                fe;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);

  //ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"Hello from process: %d\n", rank);CHKERRQ(ierr);

  ierr = processUserOptions(PETSC_COMM_WORLD, &user);CHKERRQ(ierr);
  ierr = dmMeshSetup(PETSC_COMM_WORLD, &user, &dm);CHKERRQ(ierr);
  //ierr = drawOneElem(dm,&user);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm, &user);CHKERRQ(ierr);

  ierr = DMGetApplicationContext(dm, &fe);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"fe->polydegree: %d\n", fe->polydegree);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm, &exactSol);CHKERRQ(ierr);
  //VecView(u,PETSC_VIEWER_STDOUT_WORLD);

  ierr = computeExact(dm, exactSol);
  //VecView(exactSol,PETSC_VIEWER_STDOUT_WORLD);


  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}//end of main
