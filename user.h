#include<petscdmplex.h>
/*for dm and dmplex objects*/

typedef struct ApplicationContext_private *AppCtx;

struct ApplicationContext_private{
  char          filename[PETSC_MAX_PATH_LEN]; /*optional exodusII file*/
  PetscBool	    interpolate;                  /*Generate intermediate mesh points: nodes on faces, edges in addition to vertices and element*/
  PetscInt      polydegree;          /*number of quadrature points in 1D*/
  PetscInt      dof;                          /*number of degrees of freedom at each node, e.g. For 3D elasticity: at each node we have u1, u2 and u3 */
};

PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx userOptions);
PetscErrorCode dmMeshSetup(MPI_Comm comm, AppCtx user, DM *dm);
// Use this function after calling dmMeshSetup 
PetscErrorCode drawEleme(MPI_Comm comm, DM *dm);