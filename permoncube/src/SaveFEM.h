
#ifndef SAVEFEM_H_
#define SAVEFEM_H_

#include "utility.h"
#include "Fem.h"
#include "Data.h"
#include "KSparse.h"

class CSaveFEM {
public:
  MPI_Comm comm;
private:
	CFem* fem;
	CData* data;
	CDomain* domainG;
  int MPIrank, MPIsize;

public:
	CSaveFEM(MPI_Comm comm, CFem* fem, CData* data, CDomain *domainG);
	virtual ~CSaveFEM();

private:
	void saveElements();
	void saveCoordinates();
	void saveDisplacement(double *u,double *f);
	void saveLocalStifnessMat(CKSparse * Ksparse);

public:
  void GatherDataFemToMaster();
	void saveVTK();
	void save_data();
  void saveStiffnessMatrix();
	void saveVector_int(int * x, int n, char * name);
	void saveVector_double(double * x, int n, char * name);
};

#endif /* SAVEFEM_H_ */
  
