#include "user.h"

#undef __FUNCT__
#define __FUNCT__ "processUserOptions"
PetscErrorCode processUserOptions(MPI_Comm comm, AppCtx *userOptions)
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
PetscErrorCode dmMeshSetup(MPI_Comm comm, AppCtx *user, DM *dm)
{
  PetscInt pStart, pEnd,   //all points (nodes) in mesh
           cStart,cEnd,    //cells (elements)
           vStart,vEnd,    //vertices
           dim,
	         i,j,
           numPoints,      //number of points in TransitiveClosure
           *points;        //array of points in TransitiveClosure
	const char     *filename    = user->filename;
	PetscBool      interpolate  = user->interpolate;
	PetscInt       dof          = user->dof;
	PetscErrorCode ierr;
	Vec            coords;
	DM distributedMesh = NULL;
  PetscSection   section;
  PetscInt       *conn, sz_conn;
  PetscInt       *perm_idx,         //permutation array index for local ordering
                 sz_perm_idx=0,    //size of permutation array index
                 perm27_idx[] = {22,10,19,9,1,7,21,8,20,15,3,16,5,0,6,18,4,17,24,11,23,12,2,14,25,13,26},
                 perm8_idx[] = {1,0,2,3,5,4,6,7};
  FE             fe;


	PetscFunctionBeginUser;

  ierr = PetscNew(&fe);CHKERRQ(ierr);
  fe->polydegree = user->polydegree;
  //initialize finite element space (fe)
  ierr = FESetup(fe);CHKERRQ(ierr);

  //Create a dmplex from Exodus-II file
  ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm);CHKERRQ(ierr);

  /* How to deal with boundary
  PetscInt vtx, boundary=0;
  for all vertices{                                   400, 500, 600
  ierr = DMPlexGetLabelValue(dm, "Vertex Sets",vtx, &boundary )
                                                   100 or 200, 300
  ierr = DMPlexGetLabelValue(dm, "Face Sets",fce, &boundary )
  }
  */

	/* Distribute mesh over processors */
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
  //set dof for all nodes in pStart to pEnd
  if(interpolate){  //if interpolated set section fr all nodes
    for(i=pStart; i < pEnd; i++){
      ierr = PetscSectionSetDof(section, i, dof);CHKERRQ(ierr);
    }
  }
  else{  //othersie set sectionf or vertices only
    ierr = DMPlexGetDepthStratum(*dm, 0, &vStart,&vEnd);CHKERRQ(ierr);
    for(i=pStart; i < pEnd; i++){
      if(i >= vStart && i <=vEnd){
        ierr = PetscSectionSetDof(section, i, dof);CHKERRQ(ierr);
      }
      else
      {
        ierr = PetscSectionSetDof(section, i, 0);CHKERRQ(ierr);
      }
    }
  }
  //Calculate offsets based upon the number of degrees of freedom for each point.
	ierr = PetscSectionSetUp(section);CHKERRQ(ierr);
  //Set the section for dm
	ierr = DMSetDefaultSection(*dm,section);CHKERRQ(ierr);


  if(interpolate){
    sz_perm_idx = 27;
    ierr = PetscMalloc1(sz_perm_idx,&perm_idx);CHKERRQ(ierr);
    for(i=0; i<sz_perm_idx; i++)
      perm_idx[i] = perm27_idx[i];
  }
  else{
    sz_perm_idx = 8;
    ierr = PetscMalloc1(sz_perm_idx,&perm_idx);CHKERRQ(ierr);
    for(i=0; i<sz_perm_idx; i++)
      perm_idx[i] = perm8_idx[i];
  }

  //in hierarchy of interpolated mesh topology, get all the cells (elements). (Height 0)
  ierr = DMPlexGetHeightStratum(*dm, 0, &cStart,&cEnd);CHKERRQ(ierr);
  sz_conn = sz_perm_idx*(cEnd-cStart);
  ierr = PetscMalloc1(sz_conn, &conn); CHKERRQ(ierr);
  //for each element

  for(i = cStart ; i <cEnd; i++)
  {
     points=NULL;
     //get the TransitiveClosure of each element
     ierr = DMPlexGetTransitiveClosure(*dm, i, PETSC_TRUE, &numPoints , &points);CHKERRQ(ierr);
     PetscInt k =0;
     for(j = 0; j < numPoints ; j++){
         PetscInt tmpdof;
         PetscSectionGetDof(section, points[2*j], &tmpdof);
         if (!tmpdof) continue;
         //permute the TransitiveClosure of each elemet according to perm_idx
         conn[numPoints*i+k] = points[2*perm_idx[k]];
         k=k+1;
     }
        ierr = DMPlexRestoreTransitiveClosure(*dm, i , PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
  }

    //set the connectivity (conn) to user->conn
     user->conn = conn;




    /*
     ierr = PetscPrintf(PETSC_COMM_SELF,"sz_conn %d:  \n",sz_conn);CHKERRQ(ierr);
     for(i=0; i < sz_conn; i++){
       ierr = PetscPrintf(PETSC_COMM_SELF,"conn[%d]: %d \n",i,conn[i]);CHKERRQ(ierr);
     }
    */

  //Destroy PetscSection
	PetscSectionDestroy(&section);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "drawOneElem"
PetscErrorCode drawOneElem(DM dm, AppCtx *user){

  PetscInt vStart,vEnd,    //vertices
           eStart,eEnd,    //edges
           fStart, fEnd,   //faces
           cStart,cEnd,    //cells (elements)
           rank,
           i,j,
           numPoints,      //number of points in TransitiveClosure
           *points;        //array of points in TransitiveClosure
  PetscInt const *myCone;
  PetscErrorCode ierr;
  PetscBool interpolated = user->interpolate;



  	PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);


    //in hierarchy of mesh topology, get all the vertices. (Depth 0)
    if(interpolated){
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart,&vEnd);CHKERRQ(ierr);
  	ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of vertex %d\n\n\n", rank,vEnd-vStart);CHKERRQ(ierr);
  	for(i = vStart ; i <vEnd; i++)
  	{
  		ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
  		ierr = PetscPrintf(PETSC_COMM_SELF,"process %d vertex Idx %d\n",rank, i);CHKERRQ(ierr);
  	}


    //in hierarchy of mesh topology, get all the edges. (Depth 1)
    ierr = DMPlexGetDepthStratum(dm, 1, &eStart,&eEnd);CHKERRQ(ierr);
  	ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of edges %d\n\n\n", rank,eEnd-eStart);CHKERRQ(ierr);
  	for(i = eStart ; i <eEnd; i++)
  	{
  		ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
  		ierr = PetscPrintf(PETSC_COMM_SELF,"process %d edge Idx %d : [%d,  %d]\n",rank, i, myCone[0], myCone[1]);CHKERRQ(ierr);
  	}

    //in hierarchy of mesh topology, get all the faces. (Depth 2)
    ierr = DMPlexGetDepthStratum(dm, 2, &fStart,&fEnd);CHKERRQ(ierr);
  	ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of Faces %d\n\n\n", rank,fEnd-fStart);CHKERRQ(ierr);
   	for(i = fStart ; i <fEnd; i++)
  	{
  		ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
  		ierr = PetscPrintf(PETSC_COMM_SELF,"process %d face Idx %d : [%d,  %d, %d, %d]\n",rank, i, myCone[0], myCone[1],myCone[2], myCone[3] );CHKERRQ(ierr);
  	}

    //in hierarchy of mesh topology, get all the cells (elements). (Depth 3)
    ierr = DMPlexGetDepthStratum(dm, 3, &cStart,&cEnd);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of cells %d\n\n\n", rank,cEnd-cStart);CHKERRQ(ierr);
  	for(i = cStart ; i <cEnd; i++)
  	{
  		ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
  		ierr = PetscPrintf(PETSC_COMM_SELF,"process %d cell Idx %d : [%d,  %d, %d, %d, %d, %d]\n",rank, i, myCone[0], myCone[1],myCone[2], myCone[3],myCone[4], myCone[5] );CHKERRQ(ierr);
    }

  }
  else{

    ierr = DMPlexGetDepthStratum(dm, 0, &vStart,&vEnd);CHKERRQ(ierr);
  	ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of vertex %d\n\n\n", rank,vEnd-vStart);CHKERRQ(ierr);
  	for(i = vStart ; i <vEnd; i++)
  	{
  		ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
  		ierr = PetscPrintf(PETSC_COMM_SELF,"process %d vertex Idx %d\n",rank, i);CHKERRQ(ierr);
  	}

    //in hierarchy of mesh non-interpolated topology, get all the cells (elements). (Depth 1)
    ierr = DMPlexGetDepthStratum(dm, 1, &cStart,&cEnd);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"process %d number of cells %d\n\n\n", rank,cEnd-cStart);CHKERRQ(ierr);
    for(i = cStart ; i <cEnd; i++)
    {
      ierr = DMPlexGetCone(dm, i, &myCone);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"process %d cell Idx %d : [%d,  %d, %d, %d, %d, %d, %d, %d]\n",rank, i, myCone[0], myCone[1],myCone[2], myCone[3],
                                                                                                          myCone[4], myCone[5],myCone[6], myCone[7] );CHKERRQ(ierr);
    }


  }

    for(i = cStart ; i <cEnd; i++)
    {
        points=NULL;
        //get the TransitiveClosure of each element
        ierr = DMPlexGetTransitiveClosure(dm, i, PETSC_TRUE, &numPoints , &points);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF,"\nnumber of points in TransitiveClosure:  %d\n", numPoints);CHKERRQ(ierr);
        ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,"TransitiveClosure ordering\n");CHKERRQ(ierr);
        ierr = PetscSynchronizedPrintf(PETSC_COMM_SELF,"----------------------------\n");CHKERRQ(ierr);
        for(j = 0; j < numPoints ; j++){
          ierr = PetscPrintf(PETSC_COMM_SELF,"%d--> %d\n", j, points[2*j]);CHKERRQ(ierr);
        }
          ierr = DMPlexRestoreTransitiveClosure(dm, i , PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
    }

    // NOTE 2: The code contained in the following braces {} was supposed to createa .vtu file. However, the .vtu file
    //      that it creates is faulty as the source code requires a PetscDS while we don't have it.
    // NOTE 1 : Works with 1 dof only
    // {
    //   PetscViewer view;
    //   Vec X;
    //   PetscScalar *x, *arr;
    //   PetscInt sz_X;
    //
    //   ierr = DMCreateGlobalVector(dm,&X);CHKERRQ(ierr);
    //
    //   ierr = VecGetArray(X, &x);CHKERRQ(ierr);
    //   ierr = VecGetSize(X, &sz_X); CHKERRQ(ierr);
    //   ierr = PetscMalloc1(sz_X,&arr);CHKERRQ(ierr);
    //
    //   // ierr = VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //
    //   for(i = 0 ; i<sz_X ; i++){
    //     arr[i] = i;
    //     x[user->conn[i]] += arr[i];
    //   }
    //
    //   ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
    //   ierr = PetscViewerCreate(PETSC_COMM_SELF,&view);CHKERRQ(ierr);
    //   ierr = PetscViewerSetType(view, PETSCVIEWERVTK);CHKERRQ(ierr);
    //   ierr = PetscViewerFileSetName(view, "elem.vtu");CHKERRQ(ierr);
    //   ierr = VecView(X, view);CHKERRQ(ierr);
    //   //ierr = VecView(X, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //   ierr = PetscViewerDestroy(&view);CHKERRQ(ierr);
    //   ierr = PetscFree(x);CHKERRQ(ierr);
    //   ierr = PetscFree(arr);CHKERRQ(ierr);
    // }

      PetscFunctionReturn(0);
}
