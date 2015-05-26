#ifndef DATA_H_
#define DATA_H_

#include "utility.h"
#include "KSparse.h"

class CData {
public:
  MPI_Comm comm;
	CKSparse* KSparse;
    CBSparse* B;
	CSolid45EqNumbGlobal *stif_glob_number;
	double* u;
	double* du;
	double* ddu;
	double* f;
	double* fGlobal;
	double* fE;
	double* Ku;
	double * f_subdom;
	double * f_subdom_internal;
//TODO next line changed by AM:
// check, if it does not collide
//  double * u_restrict;
 
  longInt  * indExterDOFs;
  int      n_indExterDOFs;
  

public:
	CData(MPI_Comm comm, int _i_domOnClust);
	virtual ~CData();

	void initialize(MPI_Comm comm ,int neq, int neqAll, int max_nDOFs_u_is_stored, 
                    int n_elementsSub, int n_extDOFs);
	void extDofs(CFem *fem, CDomain *domainG);
	void extDofs_l2g(int nExtDofs,std::vector<longInt> l2g);
	void copyDataToNextIter(int n_elementsSub);

};

#endif /* DATA_H_ */
