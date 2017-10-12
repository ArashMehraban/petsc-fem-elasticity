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
#include "topology.h"

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
  Topology          topo;
  PetscInt          sz_Ul,e;
  const PetscScalar *u;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);

  ierr = processUserOptions(PETSC_COMM_WORLD, &user);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_SELF,"user.ne: %d\n",user.ne);
  ierr = dmCreate(PETSC_COMM_WORLD, user, &dm);CHKERRQ(ierr);
  //ierr = drawOneElem(dm,user);CHKERRQ(ierr);

  ierr = DMGetApplicationContext(dm, &topo);CHKERRQ(ierr);

   ierr = PetscPrintf(PETSC_COMM_SELF,"topo->sz_connQ1: %d\n", topo->sz_connQ1);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF,"topo->sz_connQ2: %d\n", topo->sz_connQ2);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF,"topo->sz_perm_idx_Q1: %d\n", topo->sz_perm_idx_Q1);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_SELF,"topo->sz_perm_idx_Q2: %d\n", topo->sz_perm_idx_Q2);CHKERRQ(ierr);

  //ierr = DMCreateGlobalVector(dm, &exactSol);CHKERRQ(ierr);

  //ierr = computeExact(dm, exactSol);
  //VecView(exactSol,PETSC_VIEWER_STDOUT_WORLD);

  //ierr = DMCreateGlobalVector(dm, &res);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = VecGetLocalSize(Ul,&sz_Ul);CHKERRQ(ierr);
  //NOTE: Must populate local vector Ul from Glodbal Vector U
  // // ierr = DMGlobalToLocalBegin(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  // // ierr = DMGlobalToLocalEnd(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Ul,&u);CHKERRQ(ierr);
  //VecGetLocalSize(Ul,&sz_Ul);
  ierr = PetscPrintf(PETSC_COMM_SELF,"sz_Ul %d\n",sz_Ul);CHKERRQ(ierr);
  //VecView(Ul,PETSC_VIEWER_STDOUT_SELF);

 //test for Q1 and Q2 dmRestrictElems
  if(user.interpolate){
    PetscInt numElems;

    ierr = getNumElems(topo, &numElems);CHKERRQ(ierr);
    for(e=0; e < numElems; e+=user.ne){
      PetscScalar ue[user.dof * topo->sz_perm_idx_Q2 *  user.ne]_align;
      ierr = dmRestrictElems(dm, u, e, user.ne, Q2, user.dof ,ue);CHKERRQ(ierr);
    }

    for(e=0; e < numElems; e+=user.ne){
      PetscScalar ue[user.dof * topo->sz_perm_idx_Q1 *  user.ne]_align;
      ierr = dmRestrictElems(dm, u, e, user.ne, Q1, user.dof , ue);CHKERRQ(ierr);
    }

  }
  else
  {
    PetscInt numElems;
    ierr = getNumElems(topo, &numElems);CHKERRQ(ierr);

    for(e=0; e < numElems; e+=user.ne){
      PetscScalar ue[user.dof * topo->sz_perm_idx_Q1 *  user.ne]_align;
      ierr = dmRestrictElems(dm, u, e, user.ne, Q1, user.dof , ue);CHKERRQ(ierr);
    }
  }
  ierr = VecRestoreArrayRead(Ul,&u);CHKERRQ(ierr);


  ierr= dmFEcreate(user.dof, user.polydegree, user.addquadpts, &fe);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"fe->dof: %d\n", fe->dof);CHKERRQ(ierr);

  // ierr= DmGetCoordianteDM(dm,&dmx);CHKERRQ(ierr);
  // ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  // ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);


  //ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  //ierr = VecDestroy(&exactSol);CHKERRQ(ierr);
  ierr = VecDestroy(&Ul);CHKERRQ(ierr);



  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}//end of main
