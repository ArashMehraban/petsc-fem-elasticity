#ifndef __EXACT__H
#define __EXACT__H

#include<petscsnes.h>
#include "petscdmplex.h"
#include <petscmath.h>

PetscErrorCode givenSol(PetscReal xp[], PetscReal cp[]); //xp = u , cp = coordinate. See implemnetation
PetscErrorCode computeExact(DM dm, Vec sol);

#endif //end of __FE__H
