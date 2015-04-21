#ifndef CLUST_G_H 
#define CLUST_G_H

#include "Fem.h"
#include "Data.h"
#include "Domain_g.h"
#include "BSparseC.h"
#include "utility.h"

class CClust_g{
public:
  vector<CFem*>  fem;
  vector<CData*>  data;
  MPI_Comm comm;
  CDomain *domainG;
  int *indExterDOFsOnClust;
  int *neighbClst;
  double *u_restrict;
  CBSparseC *B1;

  int argc_l;
  char** argv_l;

  CClust_g(MPI_Comm comm,int argc, char** argv); 
  
	virtual ~CClust_g();

  void collectingExternalDOFs();
  void getGlobalIndexOfSubdomMeshGen();
  void GatherDataFemToMaster();
  void createVTK_per_cluster();
  void createVTK_per_cluster_new();
  void getMeshCoordsOfSubdom(int globoalIndSub,int *Isub, int *Jsub, int *Ksub);
  void getMeshCoordsOfClust(int globoalIndClst, int *I, int *J, int *K);



};


#endif /* CLUST_G_H*/ 
