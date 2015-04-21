#ifndef SOLVER_H_ #define SOLVER_H_

#include "utility.h"
#include "Data.h"
#include "Fem.h"
#include "KSparse.h"
#include "Solid45EqNumbGlobal.h"
#include "Cluster_g.h"
#ifndef WIN32
	#include <getopt.h>
#endif

class CSolver {
public:
	CSolver();
	virtual ~CSolver();
	static void definitionAndSolving(CData *data,	CFem *fem, CDomain *domainG,  CClust_g & clust);
	static void solver_simple_cg(int neqSub, double *u,double *f,const double &tol,const int MaxIter,CKSparse *Ksparse,CFem * fem, int verbose);
	

	static void solver_hfeti(CClust_g & clust_g);  //(CData *data, CFem*fem, CDomain *domainG);



#ifdef FLLOP_ENABLED
  static void solver_fllop(CData *data, CFem*fem, CDomain *domainG);
#endif
	static void solving(CFem *fem, CDomain *domainG, CData *data, CClust_g & clust);

};

#endif /* SOLVER_H_ */
