#include "fe.h"

#undef __FUNCT__
#define __FUNCT__ "FEbasisEval"
PetscErrorCode FEbasisEval(FE fe, PetscReal q, PetscReal B[], PetscReal D[])
{
  PetscFunctionBeginUser;
  switch(fe->polydegree){
    case 1:
      B[0] = (1-q)/2;
      B[1] = (1+q)/2;
      if(D){
        D[0] = -0.5;
        D[1] = 0.5;
      }
      break;
    case 2:
      B[0] = 0.5*(PetscSqr(q) - q); //Note: PetscSqr means Square of (NOT square root of)
      B[1] = 1 - PetscSqr(q);
      B[2] = 0.5*(PetscSqr(q) + q);
      if(D){
        D[0] = q - 0.5;
        D[1] = -2*q;
        D[2] = q + 0.5;
      }
      break;
    default: SETERRQ1(fe->comm, PETSC_ERR_SUP, "fe->polydegree %D", fe->polydegree);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FESetup"
PetscErrorCode FESetup(FE fe)
{
  PetscErrorCode ierr;
  PetscInt i,j,k,P,Q;

  PetscFunctionBeginUser;
  // Reference Element Evaluation
  P = fe->polydegree+1;
  Q = P + fe->addquadpts;  // Number of additonal quadrature points in 1D

  //Allocate space for B, D, x ,w
  //                        x: Guass points    w: Gauss weights
  ierr = PetscMalloc6(P*Q,&(fe->ref.B), P*Q,&fe->ref.D, Q*Q,&fe->ref.D_tilda, Q,&fe->ref.x, Q,&fe->ref.w, Q*Q*Q,&fe->ref.w3);CHKERRQ(ierr);
  ierr = PetscDTGaussQuadrature(Q,-1,1,fe->ref.x, fe->ref.w);CHKERRQ(ierr);
  //populate w3 which w (kron) w (kron) w
  for(i=0; i<Q; i++){
    for(j=0; j<Q; j++){
      for(k=0; k<Q; k++){
        fe->ref.w3[(i*Q+j)*Q+k] = fe->ref.w[i] * fe->ref.w[j] * fe->ref.w[k];
      }
    }
  }

  //populate B and D from Gauss quadrature points and weights;
  for(i=0; i<Q; i++){
    const PetscReal q = fe->ref.x[i];
    ierr = FEbasisEval(fe, q, &fe->ref.B[i*P], &fe->ref.D[i*P]);CHKERRQ(ierr);
  }

  if(P==2 && Q==3){

    //get the PsuedoInverse of B
    PetscReal *Bsource, *Binv;
    PetscScalar *tau, *work;
    ierr = PetscMalloc4(Q*P, &Bsource, P*Q, &Binv, Q, &tau, Q*P,&work);CHKERRQ(ierr);

    for(i=0;i<P*Q;i++)
      PetscPrintf(PETSC_COMM_SELF, "fe->ref.B[%d]: %f\n",i, fe->ref.B[i]);

    for(i=0; i<P*Q;i++)
      Bsource[i] = fe->ref.B[i];

    ierr = PetscDTPseudoInverseQR(Q, Q, P, Bsource, Binv, tau, Q, work);CHKERRQ(ierr);

    for(i=0;i<P*Q;i++)
      PetscPrintf(PETSC_COMM_SELF, "Binv[%d]:%f\n",i, Binv[i]);

    PetscScalar matBinv[P][Q];

    //On the stack, populate Matrix Binv as appose to array Binv
    for(i=0; i<P; i++){
      for(j=0; j<Q; j++){
        matBinv[i][j] = Binv[i*Q+j];
      }
    }

     PetscScalar matD[Q][P];
      //On the stack, populate Matrix D as appose to array D
      for (i=0; i<Q; i++) {
        for (j=0; j<P; j++) {
          matD[i][j] = fe->ref.D[i*P+j];
        }
      }

     // populate D_tilda = matD * matBinv  (D_tilda is in array format)
      PetscInt m=0;
      for(i=0; i<Q; i++){
       for(j=0; j<Q;j++){
         for(k=0; k<P; k++){
           fe->ref.D_tilda[m] += matD[i][k] * matBinv[k][j];
         }
         m++;
       }
      }

    //free space used in PetscDTPseudoInverseQR
    ierr = PetscFree4(Bsource, Binv, tau, work);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "dmFEcreate"
PetscErrorCode dmFEcreate(PetscInt dof, PetscInt polydegree, PetscInt addquadpts, FE *fe)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscNew(fe);CHKERRQ(ierr);
    (*fe)->dof = dof;
    (*fe)->addquadpts = addquadpts;
    (*fe)->polydegree = polydegree;
    ierr = FESetup(*fe);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
