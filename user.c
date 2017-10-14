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
  userOptions->ne = 1;              //default number of elements to extract per iteration
  userOptions->addquadpts = 0;      //default number of extra quadrature points

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

  // set the degree of polynomial for interpolation (1 less than the number of interpolating nodes in 1D)
  ierr = PetscOptionsInt("-polydegree", "Finite Element polynomial degree", "main", userOptions->polydegree, &userOptions->polydegree, NULL);CHKERRQ(ierr);

  // set the number of additional quadrature points in 1D
  ierr = PetscOptionsInt("-quad", "Finite Element polynomial degree", "main", userOptions->addquadpts, &userOptions->addquadpts, NULL);CHKERRQ(ierr);

  // set the number of quadrature points in 1D for the mesh data
  ierr = PetscOptionsInt("-dof", "number of degrees of freedom at each mesh point", "main", userOptions->dof, &userOptions->dof, NULL);CHKERRQ(ierr);

  ierr = PetscOptionsInt("-ne", "number of elements to extract per iteration", "main", userOptions->ne, &userOptions->ne, NULL);CHKERRQ(ierr);

  //End setting up user options
  PetscOptionsEnd();

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dmCreate"
PetscErrorCode dmCreate(MPI_Comm comm, AppCtx user, DM *dm)
{
  PetscInt       pStart, pEnd,   //all points (nodes) in mesh
                 cStart,cEnd,    //cells (elements)
                 numCells,       //number of cell for the current processor
                 vStart,vEnd,    //vertices
                 numPad,         //number of elements that must be padded
	               i,j,k,
                 numPoints,      //number of points in TransitiveClosure
                 *points;        //array of points in TransitiveClosure
	const char     *filename    = user.filename;
	PetscBool      interpolate  = user.interpolate;
	PetscInt       dof          = user.dof;
  PetscInt       ne           = user.ne;
	PetscErrorCode ierr;
	DM distributedMesh = NULL;
  PetscSection   section;
  PetscInt       *conn, *tmp_conn, sz_conn;
  PetscInt       *perm_idx,         //permutation array index for local ordering
                 sz_perm_idx=0,     //size of permutation array index
                 perm27_idx[] = {22,10,19,9,1,7,21,8,20,15,3,16,5,0,6,18,4,17,24,11,23,12,2,14,25,13,26},
                 perm8_idx[]  = {1,0,2,3,5,4,6,7};
  //FE             fe;
  Topology       topo;
  //Also build Q1 connectivity if Q2 is built
  PetscInt       *connQ1,
                  sz_connQ1,
                  sz_perm_idx_Q1;



	PetscFunctionBeginUser;

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
		DMDestroy(dm);
		*dm  = distributedMesh;
	}

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
  else{  //otherwise set section for vertices only
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
  numCells = cEnd-cStart;
  numPad = ne - (numCells % ne);
  if(numPad == ne)
    numPad = 0;

  // ierr = PetscPrintf(PETSC_COMM_SELF,"numCells: %d\n",numCells);CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_SELF,"numPad: %d\n",numPad);CHKERRQ(ierr);

  sz_conn = sz_perm_idx*(numCells+numPad);
  ierr = PetscCalloc1(sz_conn, &tmp_conn); CHKERRQ(ierr);

  //for each element
  for(i = cStart ; i <cEnd; i++)
  {
     points=NULL;
     //get the TransitiveClosure of each element
     ierr = DMPlexGetTransitiveClosure(*dm, i, PETSC_TRUE, &numPoints , &points);CHKERRQ(ierr);
     //ierr = PetscPrintf(PETSC_COMM_SELF,"numPpoints: %d\n",numPoints);CHKERRQ(ierr);
     PetscInt tmpOffset[27];
     PetscInt k =0;
     for(j = 0; j < numPoints ; j++){
         PetscInt tmpdof;

         PetscSectionGetDof(section, points[2*j], &tmpdof);
         if (!tmpdof)  continue;
         PetscSectionGetOffset(section, points[2*j], &tmpOffset[k++]);
       }
       for(j=0; j<k; j++){
         //permute the TransitiveClosure of each elemet according to perm_idx
         tmp_conn[sz_perm_idx*i+j] = tmpOffset[perm_idx[j]];
       }
        ierr = DMPlexRestoreTransitiveClosure(*dm, i , PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
  }

  for(i = 0 ; i <sz_conn; i++)
  {
      ierr = PetscPrintf(PETSC_COMM_SELF,"tmp_conn[%d]: %d\n",i, tmp_conn[i]);CHKERRQ(ierr);
  }



  //pad tmp_conn with the last element
  if(numPad){
    PetscInt startPad_idx, endPad_idx;

    //start to end indecies to populate with the last elements of tmp_conn
    startPad_idx = sz_perm_idx *  numCells;
    endPad_idx = sz_conn;  //NOTE: sz_conn is 1 more that the endPad_idx. So use < instead of <=

    //get the last element reordered in TransitiveClosure.
    for(i = startPad_idx; i <  endPad_idx; i++){
      tmp_conn[i] = tmp_conn[i - sz_perm_idx];
    }

  }

  //Allocate memory for conn
  ierr = PetscCalloc1(sz_conn, &conn); CHKERRQ(ierr);
  //permute tmp_conn to have it in structured by the number of elements. Store permutation in conn
  PetscInt m = 0;
  for(k = 0; k < sz_conn/(ne*sz_perm_idx);k++){
    for(j = 0; j< sz_perm_idx; j++){
      for(i = 0; i < ne; i++){
        conn[m] = tmp_conn[i*sz_perm_idx+j+k*(ne*sz_perm_idx)];
        m++;
      }
    }
  }

  //if interpolated, also build a connectivity matrix based on Q1
  if(interpolate)
  {
    PetscInt *tmp_connQ1;
    PetscInt const *myCone;

    sz_perm_idx_Q1 = 8;
    sz_connQ1 = sz_perm_idx_Q1*(numCells+numPad);

    ierr = PetscCalloc1(sz_connQ1, &tmp_connQ1); CHKERRQ(ierr);
    ierr = PetscCalloc1(sz_connQ1, &connQ1); CHKERRQ(ierr);

  	for(i = cStart ; i <cEnd; i++)
  	{
  		ierr = DMPlexGetCone(*dm, i, &myCone);CHKERRQ(ierr);
      for(j = 0; j < sz_perm_idx_Q1; j++){
        tmp_connQ1[sz_perm_idx_Q1*i+j] = myCone[perm8_idx[j]];
      }
  	}

    //permute tmp_connQ1 based on ne
    PetscInt m = 0;
    for(k = 0; k < sz_connQ1/(ne*sz_perm_idx_Q1);k++){
      for(j = 0; j< sz_perm_idx_Q1; j++){
        for(i = 0; i < ne; i++){
          connQ1[m] = tmp_connQ1[i*sz_perm_idx_Q1+j+k*(ne*sz_perm_idx_Q1)];
          m++;
        }
      }
    }
    ierr = PetscFree(tmp_connQ1);CHKERRQ(ierr);
  }

  // //setup FE
  // ierr = PetscNew(&fe);CHKERRQ(ierr);
  // fe->polydegree = user.polydegree;
  // fe->dof = user.dof;

  // //initialize finite element space (fe)
  // ierr = FESetup(fe);CHKERRQ(ierr);

  ierr = PetscNew(&topo);CHKERRQ(ierr);

  //set the connectivities to topo
  if(interpolate){
    topo->connQ2 = conn;
    topo->sz_connQ2 = sz_conn;
    topo->sz_perm_idx_Q2 = sz_perm_idx;
    topo->connQ1 = connQ1;
    topo->sz_connQ1 = sz_connQ1;
    topo->sz_perm_idx_Q1 = sz_perm_idx_Q1;
  }
  else{
    topo->connQ2 = NULL;
    topo->sz_connQ2 = 0;
    topo->sz_perm_idx_Q2 = 0;
    topo->connQ1 = conn;
    topo->sz_connQ1 = sz_conn;
    topo->sz_perm_idx_Q1 = sz_perm_idx;
  }

   // set topo to dm
   ierr = DMSetApplicationContext(*dm, topo);CHKERRQ(ierr);

  //Destroy PetscSection
	ierr = PetscSectionDestroy(&section);CHKERRQ(ierr);
  ierr = PetscFree(tmp_conn);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "drawOneElem"
PetscErrorCode drawOneElem(DM dm, AppCtx user){

  PetscInt       vStart,vEnd,    //vertices
                 eStart,eEnd,    //edges
                 fStart, fEnd,   //faces
                 cStart,cEnd,    //cells (elements)
                 rank,
                 i,j,
                 numPoints,      //number of points in TransitiveClosure
                 *points;        //array of points in TransitiveClosure
  PetscInt const *myCone;
  PetscErrorCode ierr;
  PetscBool      interpolated = user.interpolate;
  PetscScalar    *coordsArray;  //placehoder array for coordinates as a vector
  Vec            coords;        //coordiantes from DMGetCoordinate
  DM             dmc;           //For coordinates from DMGetCoordinateDM



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
        ierr = PetscPrintf(PETSC_COMM_SELF,"TransitiveClosure ordering\n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF,"----------------------------\n");CHKERRQ(ierr);
        for(j = 0; j < numPoints ; j++){
          ierr = PetscPrintf(PETSC_COMM_SELF,"%d--> %d\n", j, points[2*j]);CHKERRQ(ierr);
        }
          ierr = DMPlexRestoreTransitiveClosure(dm, i , PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
    }


    ierr = PetscPrintf(PETSC_COMM_SELF,"\n\n\n");CHKERRQ(ierr);


    //Show Coordinates
    ierr = DMGetCoordinates(dm, &coords);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(dm, &dmc);CHKERRQ(ierr);
    ierr = VecGetArray(coords, &coordsArray);CHKERRQ(ierr);

    for(i = vStart; i<vEnd; i++)
    {
      PetscScalar *cp;  //cp =[x,y,z]
      ierr = DMPlexPointLocalRead(dmc,i,coordsArray,&cp);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"coords %d:  [%f, %f,  %f]\n",i-vStart,cp[0],cp[1],cp[2]);CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(coords, &coordsArray); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
