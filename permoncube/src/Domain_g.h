#ifndef DOMAIN_G_H_
#define DOMAIN_G_H_

#include "utility.h"
#include "BoundaryCondOnFace.h"

#ifndef WIN32
	#include <getopt.h>
#endif

class CDomain {
public:
    MPI_Comm comm;
  
    //vector < int > myNeighSubDomains; // my neighboring domains 
    vector < CBoundaryCondOnFace > neuCondOnFaces; // my neighboring domains 
	  double Lx;          // length of box in x direction
	  double Ly;          // length of box in y direction
	  double Lz;          // length of box in z direction

	  int Cx;             // number of clust. in x direction
	  int Cy;             // number of clust. in y direction
	  int Cz;             // number of clust. in z direction

	  int Nx;             // number of subdom. in x direction
	  int Ny;             // number of subdom. in y direction
	  int Nz;             // number of subdom. in z direction

	  int NxClst;             // number of subdom. in x direction
	  int NyClst;             // number of subdom. in y direction
	  int NzClst;             // number of subdom. in z direction

	  int nx;             // number of elem. in x direction
	  int ny;             // number of elem. in y direction
	  int nz;             // number of elem. in z direction

	  int nxSub;             // number of elem. in x direction
	  int nySub;             // number of elem. in y direction
	  int nzSub;             // number of elem. in z direction

	  double cx;          // refinement - x
	  double cy;          // refinement - y
	  double cz;          // refinement - z

	  int nCorners_X;	  // number of corners per edge and face 
	  int nCorners_Y;
	  int nCorners_Z;
	  bool flag_DP_inner;
	  bool flag_DP_eges;

	  int n_elementsAll;  // number of all elements
	  int *n_elementsSub;  /* TODO 1d array */
	  int n_nodsAll;			// number of all coordinates
	  int *n_nodsSub;      /* TODO 1d array */
	  int n_facesAll;			// number of all edges
    int n_neighbClst;
	  int *n_facesSub;     /* TODO 1d array */
	  int neqAll;					// number of all equations
	  int *neqSub;         /* TODO 1d array */
	  int n_subdomains;		// number of subdomains
	  int n_clusters;		// number of subdomains
	  int n_subdomOnClust;		// number of subdomains
    int *n_nodOnEdg;    /* ??? */
		int *n_exterDOFs;   /* TODO 1d array */
		int n_exterDOFsClust;   /* TODO 1d array */
		int n_neighbSub;    
    int max_nDOFs_u_is_stored;

	  double fz_total;

//	  int numberOfNodesOnSubdomain;
//
    int vtk_min_rank;
    int vtk_max_rank;
    int  flag_i_callingSolver;
    bool flag_contact;
	  int flag_store_VTK;
	  bool flag_plotSol;
	  bool flag_storeCompleteMatrices;
	  bool flag_saveAllDataSendingToFortran;
	  bool flag_linear_system;
    bool flag_redund_lagr_mult;

    int verbose;
    int MPIrank;

	  double eps0;
	  double eps1;
	  int maxIter;

	  int i_load_step;
	  int i_sub_step;
	  int n_sub_steps;
	  int n_load_steps;
	  double del_deltaVec;
	  double youngsModulus;
	  double poissonsRatio;
	  double acceleration[3];
	  double density;
	  double *deltaVec;
    int *vec_neqSub;
    int *vec_neqClst;
    int *vec_n_elementsSub;
    int *vec_n_elementsClst;
    int *vec_n_nodsSub;
    int *vec_n_nodsClst;
    int *vec_numelSub;
    int *vec_numelClst;
    int *vec_globalSubNumbering;
    int *vec_globalSubNumberingMeshGen;
    int *vec_localSubNumbering;

	vector< vector < int > > lambda_map_sub; 

public:
	  CDomain(MPI_Comm comm);
	virtual ~CDomain();

public:

  void setFromOptions(int *argc, char ***args);
	void print();

};

#endif /* DOMAIN_G_H_ */
