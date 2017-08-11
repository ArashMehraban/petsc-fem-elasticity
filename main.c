static char help[] = "fem-elasticity";
//#include<petscds.h>
#include<petscdt.h>
/* includes the header for PetscDTGaussQuadrature
 that calculates the quadrature points and weights*/
#include<petscdmplex.h>
/*for dm and dmplex objects*/
#include<petscsnes.h>
/*petscsnes.h also includes the following headers:
  1) petscsys.h
  2) petscmat.h
  3) petscis.h
  4) petscviewer.h
  5) petscksp.h
  6) petscpc.h  
*/
#include"basis.h"

typedef struct{
  char          filename[PETSC_MAX_PATH_LEN]; /*optional exodusII file*/
  PetscBool	    interpolate;                  /*Generate intermediate mesh elements*/
  PetscInt      num_quadr_pts_in_1d;          /*number of quadrature points in 1D*/
  PetscInt      dof;                          /*number of degrees of freedom at each node, e.g. For 3D elasticity: at each node we have u1, u2 and u3 */

}AppCtx;


//PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
//PetscErrorCode FormJacobian(SNES,Vec, Mat,Mat, void*);
//static PetscErrorCode kron(PetscScalar *A,PetscInt nrA,PetscInt ncA, PetscScalar *B, PetscInt nrB, PetscInt ncB, PetscScalar *C, PetscInt *nrC, PetscInt *ncC);
/*static PetscErrorCode allocateBasis(PetscInt num_quadr_pts_in_1d);
static PetscErrorCode destroyBasis();*/
//static PetscErrorCode allocateTensor(PetscInt num_quadr_pts_in_1d, PetscInt dim,PetscScalar *B, PetscScalar *D, PetscScalar *W);
//static PetscErrorCode destroyTensor(PetscScalar *B, PetscScalar *D, PetscScalar *W);
//static PetscErrorCode tensor(PetscInt flag);
static PetscErrorCode getShape(AppCtx *user, PetscInt dim);
static PetscErrorCode printLocalMat(PetscScalar *Mat,PetscInt nrow, PetscInt ncol);
static PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx *userOptions);
static PetscErrorCode getMesh(MPI_Comm comm, AppCtx *user, DM *dm);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  PetscMPIInt       rank;
  PetscErrorCode    ierr;
  DM                dm;
  AppCtx	        user; /*user-defined work context*/
  PetscInt          dim;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);



  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Hello from process: %d\n", rank);CHKERRQ(ierr);
 
  ierr = processUserOptions(PETSC_COMM_WORLD, &user);CHKERRQ(ierr);
  ierr = getMesh(PETSC_COMM_WORLD, &user, &dm);CHKERRQ(ierr);

  DMView(dm,PETSC_VIEWER_STDOUT_WORLD);

 
  PetscScalar myMat[][2] ={{1,2},{3,4},{5,6}};
  printLocalMat((PetscScalar *)myMat,3,2); 
  
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_SELF,"Mesh dimension is %d\n\n\n", dim);

  ierr = getShape(&user,dim);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
 
return 0;
}//end of main

// implementations (Move them to different file later)
#undef __FUNCT__
#define __FUNCT__ "processUserOptions"
static PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx *userOptions)
{

  PetscBool                  flg;
  PetscErrorCode			 ierr;
  
  PetscFunctionBeginUser;
  userOptions->interpolate = PETSC_FALSE;
  userOptions->filename[0] = '\0';
  userOptions->num_quadr_pts_in_1d = 0;
  userOptions->dof = 0;
  
  //start setting user options from command line
  ierr = PetscOptionsBegin(comm, "", "FEM elasticity options such as mesh filename, number of quadrature points in 1D, etc. ", "DMPLEX");CHKERRQ(ierr);


  //set the option (-f) to get a filename from user: mpiexec -n 1 ./main  -f filename
  ierr = PetscOptionsString("-f", "Exodus.II filename to read", "main.c", userOptions->filename, userOptions->filename, sizeof(userOptions->filename), &flg);CHKERRQ(ierr);
  #if !defined(PETSC_HAVE_EXODUSII)
    if(flg)  SETERRQ(comm, PETSC_ERR_ARG_WRONG, "This option requires ExodusII support. Reconfigure your Arch with --download-exodusii");
  #endif

  /*set the interpolate option to whether generate intermediate mesh elements. if interpolate is PETSC_TRUE and you have 2D mesh,
  with elements (Cells) and nodes only, Petsc will generate faces for you as well.*/
    ierr = PetscOptionsBool("-interpolate", "Generate intermediate mesh elements", "main", userOptions->interpolate, &userOptions->interpolate, NULL);  CHKERRQ(ierr);

  // set the number of quadrature points in 1D for the mesh data
  ierr = PetscOptionsInt("-nq1d", "number of quadrature points in 1D", "main", userOptions->num_quadr_pts_in_1d, &userOptions->num_quadr_pts_in_1d, NULL);CHKERRQ(ierr);

  // set the number of quadrature points in 1D for the mesh data
  ierr = PetscOptionsInt("-dof", "number of degrees of freedom at each mesh point", "main", userOptions->dof, &userOptions->dof, NULL);CHKERRQ(ierr);
   
  //End setting up user options
  PetscOptionsEnd();

  PetscFunctionReturn(0);
}
static PetscErrorCode printLocalMat(PetscScalar *Mat, PetscInt nrow, PetscInt ncol)
{
  PetscInt       i,j;
  PetscErrorCode ierr;
  
  for(i = 0; i<nrow; i++)
  {
     for(j = 0; j<ncol; j++)
        ierr = PetscPrintf(PETSC_COMM_SELF,"%g ", Mat[i*ncol+j]);CHKERRQ(ierr);
     ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
  }
 return 0;
}

/*#undef __FUNCT__
#define __FUNCT__ "allocateBasis"
static PetscErrorCode allocateBasis(PetscInt num_quadr_pts_in_1d)
{
	PetscErrorCode ierr;
	PetscScalar *bhat, *dhat;
	PetscInt size = num_quadr_pts_in_1d*num_quadr_pts_in_1d;  //for bhat, dhat

	PetscFunctionBeginUser;
	//Allocate space for bhat and dhat
	  ierr = PetscMalloc1(size,&bhat);CHKERRQ(ierr);
	  ierr = PetscMalloc1(size,&dhat);CHKERRQ(ierr);
	PetscFunctionReturn(0);

}*/

//#undef __FUNCT__
//#define __FUNCT__ "allocateTensor"
/*static PetscErrorCode allocateTensor(PetscInt num_quadr_pts_in_1d, PetscInt dim,PetscScalar *B, PetscScalar *D, PetscScalar *W)
{
	PetscErrorCode ierr;
	PetscInt size = num_quadr_pts_in_1d*num_quadr_pts_in_1d;

	PetscFunctionBeginUser;
      //Allocate space for tensor B and D
	  ierr = PetscMalloc1(size*size,&B);CHKERRQ(ierr);
	  ierr = PetscMalloc1(size*size*dim,&D);CHKERRQ(ierr);
	  ierr = PetscMalloc1(size,&W);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
*/
//#undef __FUNCT__
//#define __FUNCT__ "destroyTensor"
/*static PetscErrorCode destroyTensor(PetscScalar *B,PetscScalar *D,PetscScalar *W)
{
	PetscFunctionBeginUser;
	PetscFree(B);
	PetscFree(D);
	PetscFree(W);
	PetscFunctionReturn(0);
}
*/
/*
#undef __FUNCT__
#define __FUNCT__ "tensor"
static PetscErrorCode tensor(PetscInt flag)
{
	PetscFunctionBeginUser;

	PetscFunctionReturn(0);


}
*/
#undef __FUNCT__
#define __FUNCT__ "getShape"
static PetscErrorCode getShape(AppCtx *user, PetscInt dim)
{
  PetscScalar *x, *w, *bhat, *dhat;
  PetscErrorCode ierr;
  PetscInt num_quadr_pts_in_1d = user->num_quadr_pts_in_1d;
  PetscInt size = num_quadr_pts_in_1d*num_quadr_pts_in_1d ,i,j;  //for bhat, dhat

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(num_quadr_pts_in_1d,&x);CHKERRQ(ierr);
  ierr = PetscMalloc1(num_quadr_pts_in_1d,&w);CHKERRQ(ierr);
 
  ierr = PetscDTGaussQuadrature(num_quadr_pts_in_1d,-1,1, x, w); CHKERRQ(ierr);
 
  //delete this----section--> 
  PetscScalarView(num_quadr_pts_in_1d,x,PETSC_VIEWER_STDOUT_SELF); 
  PetscScalarView(num_quadr_pts_in_1d,w,PETSC_VIEWER_STDOUT_SELF); 
  //delete this----section-->

  //Allocate space for bhat and dhat
  ierr = PetscMalloc1(size,&bhat);CHKERRQ(ierr);
  ierr = PetscMalloc1(size,&dhat);CHKERRQ(ierr);

  //Evaluate basis (shape) functions at quadrature points.
  if(num_quadr_pts_in_1d == 2)
  {
	  j =0;
	  for (i=0; i < size; i = i+2)
	  {
		  bhat[i]   = (1 - x[j])/2;
		  bhat[i+1] = (1 + x[j])/2;
		  dhat[i]   = -0.5;
		  dhat[i+1] = 0.5;
		  j = j+1;
	  }
  }
  if(num_quadr_pts_in_1d == 3)
  {
	  j =0;
	  	  for (i=0; i < size; i = i+3)
	  	  {
	  		  bhat[i]   = (x[j]*x[j] - x[j])/2;
	  		  bhat[i+1] = 1 - x[j]*x[j];
	  		  bhat[i+2] = (x[j]*x[j] +x[j])/2;
	  		  dhat[i]   = x[j] - 0.5;
	  		  dhat[i+1] = -2*x[j];
	  		  dhat[i+2] = x[j] + 0.5;
	  		  j = j+1;
	  	  }

  }


  //delete this----section-->
    PetscScalarView(size,bhat,PETSC_VIEWER_STDOUT_SELF);
    PetscScalarView(size,dhat,PETSC_VIEWER_STDOUT_SELF);
    //delete this----section-->


  // free local arrays and matrices to getShape function
  ierr = PetscFree(x);CHKERRQ(ierr);
  ierr = PetscFree(w);CHKERRQ(ierr);
  ierr = PetscFree(bhat);CHKERRQ(ierr);
  ierr = PetscFree(dhat);CHKERRQ(ierr);

  PetscFunctionReturn(0);

 //return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getMesh"
static PetscErrorCode getMesh(MPI_Comm comm, AppCtx *user, DM *dm)
{
	PetscInt       rank,i,cStart,cEnd,numCells, dim;
	const char     *filename    = user->filename;
	PetscBool      interpolate  = user->interpolate;
	PetscInt       dof          = user->dof;
	PetscErrorCode ierr;
	Vec            coords;
	PetscInt const *myCone;
	PetscSection   section;
	Vec            exact_u;
	PetscScalar    init_val = 0;   /*initial value for exact_u global vector, set to zero initially*/
	PetscSection s;
	PetscInt vStart,vEnd;
	PetscInt offset;
	PetscInt pstart, pend;
	DM distributedMesh = NULL;


	PetscFunctionBeginUser;
	ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm);CHKERRQ(ierr);
	//ierr = DMGetCoordinates(*dm, &coords); CHKERRQ(ierr);
	//VecView(coords,PETSC_VIEWER_STDOUT_WORLD);

	/* Distribute mesh over processes */
	DMPlexDistribute(*dm, 0, NULL, &distributedMesh);
	if (distributedMesh) {
		DMDestroy(dm);
		*dm  = distributedMesh;
	}


	// initial petsc section here
	PetscSectionCreate(comm,&s);

	PetscSectionSetNumFields(s, 1);
	DMPlexGetChart(*dm, &pstart, &pend);
	PetscSectionSetChart(s,pstart, pend);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"process %d Section Chart %d - %d\n", rank, pstart, pend);
	DMPlexGetDepthStratum(*dm, 0, &vStart,&vEnd);
	for(i=vStart; i< vEnd; i++ )
		PetscSectionSetDof(s,i,2);
	/*For sides and element nodes
	//Declare eStart, eEnd, cStart, cEnd above
	DMPlexGetDepthStratum(*dm, 1, &eStart,&eEnd);
	for(i=eStart; i< eEnd; i++ )
		PetscSectionSetDof(s,i,2);
	DMPlexGetDepthStratum(*dm, 2, &cStart,&cEnd);
	for(i=cStart; i< cEnd; i++ )
		PetscSectionSetDof(s,i,2);
        */

	PetscSectionSetUp(s);

	for(i=pstart; i< pend; i++ )
	{
		//PetscSectionGetOffset(s,i,&offset);
		//PetscPrintf(PETSC_COMM_SELF,"mesh point %d with offset %d\n", i, offset);
	}


	DMSetDefaultSection(*dm,s);

	ierr = DMCreateGlobalVector(*dm, &exact_u);CHKERRQ(ierr);
	//VecView(exact_u,PETSC_VIEWER_STDOUT_WORLD );


	ierr = DMGetDimension(*dm, &dim);
	PetscPrintf(PETSC_COMM_WORLD,"dimenstion from dm is *********  %d\n\n\n", dim);





	DMPlexGetHeightStratum(*dm, 0, &cStart, &cEnd);
	numCells = cEnd - cStart;
	DMPlexGetChart(*dm, &pstart, &pend);
	ierr = DMPlexGetCone(*dm, pstart, &myCone);
	//ierr = DMGetTransitiveClosure(*dm, &u); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"process %d number of cells %d\n\n\n", rank,numCells);
	//DMGetDefaultSection(*dm, &section);
	//PetscSectionView(section,PETSC_VIEWER_STDOUT_WORLD);

	//ierr = DMGetDimension(*dm,&dim);CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"Mesh dimension is %d\n\n\n", dim);

	for(i = cStart ; i <cEnd; i++)
	{ierr = DMPlexGetCone(*dm, i, &myCone);
	   PetscPrintf(PETSC_COMM_SELF,"process %d cell %d = [%d %d %d %d]\n",rank, i, myCone[0],myCone[1],myCone[2],myCone[3]);

	}

	PetscSectionDestroy(&s);

	PetscFunctionReturn(0);

}



/*static PetscErrorCode kron(PetscScalar *A,PetscInt nrA,PetscInt ncA, PetscScalar *B, PetscInt nrB, PetscInt ncB, PetscScalar *C, PetscInt *nrC, PetscInt *ncC)
{
  
  PetscInt r = nrA*nrB;
  PetscInt c = ncA*ncB;
  PetscErrorCode ierr;

  ierr = PetscCalloc1(r*c, &C);CHKERRQ(ierr); 
  
 return 0;
}*/
