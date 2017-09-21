#include "exact.h"

#undef __FUNCT__
#define __FUNCT__ "givenSol"
PetscErrorCode givenSol(PetscReal xp[], PetscReal cp[])
{
  //xp =[u1,u2,u3] and cp =[x,y,z]
  // x_coord = cp[0]
  // y_coord = cp[1]
  // z_coord = cp[2]
  // u1 = xp[0]
  // u2 = xp[1]
  // u3 = xp[2]

  //Assume we have the following exact solutions:
  /* u1 = exp(2*x)*sin(3*y)*cos(4*z)
     u2 = exp(3*y)*sin(4*z)*cos(2*x)
     u3 = exp(4*z)*sin(2*x)*cos(3*y)
     */
     PetscFunctionBeginUser;
  /*
     xp[0] = PetscExpReal(2*cp[0])*PetscSinReal(3*cp[1])*PetscCosReal(4*cp[2]);
     xp[1] = PetscExpReal(3*cp[1])*PetscSinReal(4*cp[2])*PetscCosReal(2*cp[0]);
     xp[2] = PetscExpReal(4*cp[2])*PetscSinReal(2*cp[0])*PetscCosReal(3*cp[1]); */

     xp[0] = cp[0]*cp[1]*cp[2];
     xp[1] = cp[0]*cp[1]*cp[2];
     xp[2] = cp[0]*cp[1]*cp[2];

     PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "computeExact"
PetscErrorCode computeExact(DM dm, Vec sol)
{
  PetscReal *solArray, *coordsArray;
  PetscInt  //pStart, pEnd,
            i, vStart,vEnd;
  PetscErrorCode ierr;
  Vec coords;

  PetscFunctionBeginUser;

  ierr = VecSet(sol,0);CHKERRQ(ierr);
  ierr = DMGetCoordinates(dm, &coords);CHKERRQ(ierr);
  //ierr = VecView(coords, PETSC_VIEWER_STDOUT_WORLD 	);CHKERRQ(ierr);
  ierr = VecGetArray(sol, &solArray);CHKERRQ(ierr);
  ierr = VecGetArray(coords, &coordsArray);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart,&vEnd);CHKERRQ(ierr);

  for(i = vStart; i<vEnd; i++)
  {
    PetscReal *xp, *cp;
    ierr = DMPlexPointLocalRef(dm,i,solArray,&xp);CHKERRQ(ierr);
    ierr = DMPlexPointLocalRead(dm,i,coordsArray,&cp);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF,"coords %d:  [%f, %f,  %f]\n",i,cp[0],cp[1],cp[2]);CHKERRQ(ierr);
    ierr = givenSol(xp, cp);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(sol, &solArray); CHKERRQ(ierr);
  ierr = VecRestoreArray(coords, &coordsArray); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
