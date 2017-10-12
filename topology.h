#ifndef __TOPOLOGY__H
#define __TOPOLOGY__H

#include<petscdmplex.h>
#include<petscdt.h>

 typedef enum {Q1, Q2} RestricMode;

 typedef struct Topology_private *Topology;

 struct Topology_private
{
  MPI_Comm    comm;
  PetscInt    *connQ1,
              *connQ2,
              sz_connQ1,
              sz_connQ2,
              sz_perm_idx_Q1,
              sz_perm_idx_Q2;
};

PetscErrorCode dmRestrictElems(DM dm, const PetscScalar *u, PetscInt elem, PetscInt ne, RestricMode rmode, PetscInt dof, PetscScalar *y);
PetscErrorCode getNumElems(Topology topo, PetscInt *numElems);

#endif //end of __TOPOLOGY__H
