#include "topology.h"

#undef __FUNCT__
#define __FUNCT__ "dmRestrictElems"
PetscErrorCode dmRestrictElems(DM dm, const PetscScalar *u, PetscInt elem, PetscInt ne, RestricMode rmode, PetscInt dof, PetscScalar *y)
{                                                      //u is a vecArray (read-only) passed to this function
  PetscErrorCode ierr;
  Topology       topo;
  PetscInt       i, d, //for loop counters
                 *connPtr = 0,
                 elemStart, elemEnd,
                 sz_perm_idx;

  PetscFunctionBeginUser;

  ierr = DMGetApplicationContext(dm,&topo);CHKERRQ(ierr);

  if(rmode == Q1)
  {
    elemStart = topo->sz_perm_idx_Q1*elem;
    //ierr = PetscPrintf(PETSC_COMM_SELF,"elemStart: %d\n", elemStart);CHKERRQ(ierr);
    elemEnd = topo->sz_perm_idx_Q1*(elem+ne);
    connPtr = &topo->connQ1[0];
    sz_perm_idx = topo->sz_perm_idx_Q1;
  }
  else if (rmode == Q2){
    elemStart = topo->sz_perm_idx_Q2*elem;
    elemEnd = topo->sz_perm_idx_Q2*(elem+ne);
    connPtr = &topo->connQ2[0];
    sz_perm_idx = topo->sz_perm_idx_Q2;
  }
  else
  {
    SETERRQ1(topo->comm, PETSC_ERR_SUP, "RestricMode %D", rmode);
  }

  for(i = elemStart; i<elemEnd; i++){
    const PetscScalar *u_dof;

    u_dof = &u[dof*connPtr[i]];
    for(d = 0 ; d < dof; d++){
      y[d*(ne*sz_perm_idx)+i-elemStart] = u_dof[d];
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getNumElems"
PetscErrorCode getNumElems(Topology topo, PetscInt *numElems)
{
  PetscFunctionBeginUser;

  *numElems = topo->sz_connQ1 / topo->sz_perm_idx_Q1;

  PetscFunctionReturn(0);
}
