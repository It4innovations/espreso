#include "RSparse.h"

CRSparse::~CRSparse() {
  if (Rfull) delete [] Rfull;
}

CRSparse::CRSparse(CFem *fem, int neqSub) {
// TODO get kernel from   factorization
	int tmp_k;
	int inod;
	double *R = new double [6*neqSub];
  CCoordinate *coordinate = fem->mesh.coordinateSub;
	// zeroing - change it!!!
  n_row = neqSub;
  n_col = 6;
	for (int i = 0;i<6*neqSub;i++){
		R[i]=0.0;
	}
	//
	for (int i = 0;i<neqSub;i++){
		tmp_k = i%3;//fem->mesh.l2g[i]%3;
		//inod = fem->mesh.l2g[i]/3;
		inod = i/3;
		//
		if (tmp_k==0){
			R[i+0*neqSub] = 1.0;
			R[i+3*neqSub] = -coordinate[inod].y;
			R[i+4*neqSub] = -coordinate[inod].z;
		}
		else if (tmp_k==1){
			R[i+neqSub] = 1.0;
			R[i+3*neqSub] =  coordinate[inod].x;
			R[i+5*neqSub] = -coordinate[inod].z;
		}
		else if (tmp_k==2){
			R[i+2*neqSub] = 1.0;
			R[i+4*neqSub] =  coordinate[inod].x;
			R[i+5*neqSub] =  coordinate[inod].y;
		}
	}
  Rfull = R;
}
