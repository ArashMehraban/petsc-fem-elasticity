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
  PetscErrorCode    ierr;
  DM                dm;
  AppCtx	          user; /*user-defined work context*/
  Vec               exactSol, Ul, dmx ,X; //res
  FE                fe;
  PetscInt          sz_Ul,e;
  const PetscScalar *u;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);

  ierr = processUserOptions(PETSC_COMM_WORLD, &user);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_SELF,"user.ne: %d\n",user.ne);
  ierr = dmCreate(PETSC_COMM_WORLD, user, &dm);CHKERRQ(ierr);
  //ierr = drawOneElem(dm,user);CHKERRQ(ierr);

   ierr = DMGetApplicationContext(dm, &fe);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF,"fe->polydegree: %d\n", fe->polydegree);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF,"fe->dof: %d\n", fe->dof);CHKERRQ(ierr);

  //ierr = DMCreateGlobalVector(dm, &exactSol);CHKERRQ(ierr);

  //ierr = computeExact(dm, exactSol);
  //VecView(exactSol,PETSC_VIEWER_STDOUT_WORLD);

  //ierr = DMCreateGlobalVector(dm, &res);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = VecGetLocalSize(Ul,&sz_Ul);CHKERRQ(ierr);
  //NOTE: Must populate local vector Ul from Glodbal Vector U
  // ierr = DMGlobalToLocalBegin(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  // ierr = DMGlobalToLocalEnd(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Ul,&u);CHKERRQ(ierr);
  //VecGetLocalSize(Ul,&sz_Ul);
  ierr = PetscPrintf(PETSC_COMM_SELF,"sz_Ul %d\n",sz_Ul);CHKERRQ(ierr);
  //VecView(Ul,PETSC_VIEWER_STDOUT_SELF);
  for(e=0; e < (fe->sz_conn/ fe->sz_perm_idx); e+=user.ne){
    PetscScalar ue[fe->dof * fe->sz_perm_idx * user.ne]_align;
    ierr = dmExtractElems(dm, u, e, user.ne, ue);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(Ul,&u);CHKERRQ(ierr);

  // ierr= DmGetCoordianteDM(dm,&dmx);CHKERRQ(ierr);
  // ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  //ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);


  //ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = VecDestroy(&exactSol);CHKERRQ(ierr);
  ierr = VecDestroy(&Ul);CHKERRQ(ierr);



  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}//end of main
