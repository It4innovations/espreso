
#ifndef FEM_H_
#define FEM_H_

#include "utility.h"
#include "Mesh.h"
#include "Domain_g.h"
#include "BoundaryCondition.h"
#include "KSparse.h"
#include "BSparse.h"
#include "RSparse.h"
#include "LinearAlgebra.h"
#include "CustomData.h"

#include "espreso.h"

class CKSparse;
class CBSparse;
class CRSparse;
class CCustomData;

class CFem {
public:
  CMesh mesh;
  //CDomain *domainG;
  CBoundaryCondition bound_cond;
  MPI_Comm comm;
  int i_domOnClust;
//  int i_domGlobal;
//	CFem();
	CFem(const CFem& other);
  CFem & operator=(const CFem &fem){printf("in Fem 29 +++++++++\n");};
	CFem(MPI_Comm comm, int _i_domOnClust);

	virtual ~CFem();

	void element_generator3d(Mesh &mesh);
	void coordinate_generator3d(Coordinates &coordinates);
	void mesh_generator3d(CDomain *domainG);
//	void dataDirBC();
	void dataDirBCSub(int *i_face,int n_facesWithDirBC, 
						int n_facesSub, int *dir_xyz);

//	void dataConBC();
  void dataConBCSub(CDomain *domainG, int *i_face,int n_facesWithConBC, 
    				int n_facesSub, int *dir_xyz);
	void applicationDirBCtoAxb(CKSparse *Ksparse, double *f,CDomain *domainG);
	void dataNeuBC(double * f, int n, double fz_total,int i_face, int n_faces);
	void exportRHS(double *f);
	void statisticPrintToDisp(CKSparse *Ksparse, CDomain *domainG, double *u,
	        double * f, double * Ku);
  void initialize(MPI_Comm comm, int _i_domOnClust);
//	void GetMemoryStat();
  void get_DP_DOFs(CDomain *domainG);

private:   
	static double relativeForceEquilib(double *Ku, double *f, int indBC[], int nDir, int ndof);
};

#endif /* FEM_H_ */
