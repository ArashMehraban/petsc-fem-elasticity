#include "user.h"

#undef __FUNCT__
#define __FUNCT__ "processUserOptions"
PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx userOptions)
{
  PetscErrorCode	ierr;
  PetscBool       flg;

  PetscFunctionBeginUser;
  userOptions->interpolate = PETSC_FALSE;
  userOptions->filename[0] = '\0';  // mesh file to read
  userOptions->polydegree = 0;      //Finite Element polynomial degree
  userOptions->dof = 0;             //degrees of freedom (dof) per vertex

  //start setting user options from command line
  ierr = PetscOptionsBegin(comm, "", "FEM elasticity options such as mesh filename, polynomial degree,  etc. ", "DMPLEX");CHKERRQ(ierr);

  //set the option (-f) to get a filename from user: mpiexec -n 1 ./main  -f filename
  ierr = PetscOptionsString("-f", "Exodus.II filename to read", "main.c", userOptions->filename, userOptions->filename, sizeof(userOptions->filename), &flg);CHKERRQ(ierr);
  #if !defined(PETSC_HAVE_EXODUSII)
    if(flg)  SETERRQ(comm, PETSC_ERR_ARG_WRONG, "This option requires ExodusII support. Reconfigure your Arch with --download-exodusii");
  #endif

  /*set the interpolate option to generate intermediate mesh elements. Example: if interpolate is PETSC_TRUE and you have a 3D mesh, then Petsc will generate
  sides and faces for you, in addition to what is read (cells and vertices).*/
    ierr = PetscOptionsBool("-interpolate", "Generate intermediate mesh elements", "main", userOptions->interpolate, &userOptions->interpolate, NULL);  CHKERRQ(ierr);

  // set the number of quadrature points in 1D for the mesh data
  ierr = PetscOptionsInt("-polydegree", "Finite Element polynomial degree", "main", userOptions->polydegree, &userOptions->polydegree, NULL);CHKERRQ(ierr);

  // set the number of quadrature points in 1D for the mesh data
  ierr = PetscOptionsInt("-dof", "number of degrees of freedom at each mesh point", "main", userOptions->dof, &userOptions->dof, NULL);CHKERRQ(ierr);

  //End setting up user options
  PetscOptionsEnd();

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dmMeshSetup"
PetscErrorCode dmMeshSetup(MPI_Comm comm, AppCtx user, DM *dm)
{
  PetscInt pStart, pEnd,   //all points (nodes) in mesh
           vStart,vEnd,    //vertices
           eStart,eEnd,    //edges
           fStart, fEnd,   //faces
           cStart,cEnd,    //cells (elements)
	         i,
           //rank,           // <------------Delete this--------------->
           dim;            //dimension
	const char     *filename    = user->filename;
	PetscBool      interpolate  = user->interpolate;
	PetscInt       dof          = user->dof;
	PetscErrorCode ierr;
	Vec            coords;
	//PetscInt const *myCone;
	DM distributedMesh = NULL;
  PetscSection section;
  PetscInt *conn;


	PetscFunctionBeginUser;

  //Create a dmplex from Exodus-II file
  ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm);CHKERRQ(ierr);

	/* Distribute mesh over processes */
	DMPlexDistribute(*dm, 0, NULL, &distributedMesh);
	if (distributedMesh) {
		PetscPrintf(PETSC_COMM_WORLD,"Mesh is distributed\n");
		DMDestroy(dm);
		*dm  = distributedMesh;
	}

  //Set the Coordiantes for dmplex
  ierr = DMGetCoordinates(*dm, &coords);CHKERRQ(ierr);
  //Set the dimension for dmplex
  ierr = DMGetDimension(*dm, &dim);CHKERRQ(ierr);


	//Initialize petsc section
  ierr = PetscSectionCreate(comm,&section);CHKERRQ(ierr);
  //Set the number of fields from userCtx
	ierr = PetscSectionSetNumFields(section, 1);CHKERRQ(ierr);
  //get the chart from dmplex
	ierr = DMPlexGetChart(*dm, &pStart, &pEnd);CHKERRQ(ierr);
  //set the chart for each processor from pStart to pEnd
	ierr = PetscSectionSetChart(section,pStart, pEnd);CHKERRQ(ierr);

  /*  <------------Delete this Section-------------->
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_SELF,"process %d Section Chart %d - %d\n", rank, pStart, pEnd);
  */ // <------------Delete this Section-------------->

  //in hierarchy of mesh topology, get all the vertices. (Depth 0)
	ierr = DMPlexGetDepthStratum(*dm, 0, &vStart,&vEnd);CHKERRQ(ierr);
	for(i=vStart; i< vEnd; i++ )
		ierr = PetscSectionSetDof(section,i,dof);CHKERRQ(ierr);
       /*
	PetscPrintf(PETSC_COMM_SELF,"process %d number of vertex %d\n\n\n", rank,vEnd-vStart);
	for(i = vStart ; i <vEnd; i++)
	{
		ierr = DMPlexGetCone(*dm, i, &myCone);
		PetscPrintf(PETSC_COMM_SELF,"process %d vertex Idx %d\n",rank, i);
	}
       */

  //in hierarchy of mesh topology, get all the edges. (Depth 1)
	ierr = DMPlexGetDepthStratum(*dm, 1, &eStart,&eEnd);CHKERRQ(ierr);
	for(i=eStart; i< eEnd; i++ )
		ierr = PetscSectionSetDof(section,i,dof);CHKERRQ(ierr);
       /*
	PetscPrintf(PETSC_COMM_SELF,"process %d number of edges %d\n\n\n", rank,eEnd-eStart);
	for(i = eStart ; i <eEnd; i++)
	{
		ierr = DMPlexGetCone(*dm, i, &myCone);
		PetscPrintf(PETSC_COMM_SELF,"process %d edge Idx %d : [%d,  %d]\n",rank, i, myCone[0], myCone[1]);
	}
       */

  //in hierarchy of mesh topology, get all the faces. (Depth 2)
	ierr = DMPlexGetDepthStratum(*dm, 2, &fStart,&fEnd);CHKERRQ(ierr);
	for(i=fStart; i< fEnd; i++ )
		ierr = PetscSectionSetDof(section,i,dof);CHKERRQ(ierr);
       /*
	PetscPrintf(PETSC_COMM_SELF,"process %d number of Faces %d\n\n\n", rank,fEnd-fStart);
 	for(i = fStart ; i <fEnd; i++)
	{
		ierr = DMPlexGetCone(*dm, i, &myCone);
		PetscPrintf(PETSC_COMM_SELF,"process %d face Idx %d : [%d,  %d, %d, %d]\n",rank, i, myCone[0], myCone[1],myCone[2], myCone[3] );
	}
       */

  //in hierarchy of mesh topology, get all the cells (elements). (Depth 3)
	ierr = DMPlexGetDepthStratum(*dm, 3, &cStart,&cEnd);CHKERRQ(ierr);
	for(i=cStart; i< cEnd; i++ )
		ierr = PetscSectionSetDof(section,i,dof);CHKERRQ(ierr);
       /*
	PetscPrintf(PETSC_COMM_SELF,"process %d number of cells %d\n\n\n", rank,cEnd-cStart);
	for(i = cStart ; i <cEnd; i++)
	{
		ierr = DMPlexGetCone(*dm, i, &myCone);
		PetscPrintf(PETSC_COMM_SELF,"process %d cell Idx %d : [%d,  %d, %d, %d, %d, %d]\n",rank, i, myCone[0], myCone[1],myCone[2], myCone[3],myCone[4], myCone[5] );
	}
       */

	ierr = PetscSectionSetUp(section);CHKERRQ(ierr);

	ierr = DMSetDefaultSection(*dm,section);CHKERRQ(ierr);



  //for each element
  for(i = cStart ; i <cEnd; i++)
	{
    PetscInt j,
    numPoints,      //number of points in TransitiveClosure
    *points;        //array of points in TransitiveClosure

	  points=NULL;
    //get the TransitiveClosure of each element
    ierr = DMPlexGetTransitiveClosure(*dm, i, PETSC_TRUE, &numPoints , &points);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_SELF,"numPoints  %d\n\n", numPoints);CHKERRQ(ierr);
	  for(j = 0; j < numPoints ; j++)
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d--> %d\n", j, points[2*j]);CHKERRQ(ierr);
	    ierr = DMPlexRestoreTransitiveClosure(*dm, i , PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
   }




  // <-------------------Question------------->
  //Should this be destryoed here?
	PetscSectionDestroy(&section);

	PetscFunctionReturn(0);
}
