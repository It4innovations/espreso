
#include "Domain_g.h"

bool gv_flag_linear_system;

CDomain::CDomain(MPI_Comm comm):
comm(comm),

Lx(0), /*number of box length in z direction*/
Ly(0), /*number of box length in z direction*/
Lz(0), /*number of box length in z direction*/

Cx(0), /*number of clust. in x direction*/
Cy(0), /*number of clust. in y direction*/
Cz(0), /*number of clust. in z direction*/

Nx(0), /*number of subdom. in x direction*/
Ny(0), /*number of subdom. in y direction*/
Nz(0), /*number of subdom. in z direction*/

NxClst(0), /*number of subdom. in x direction*/
NyClst(0), /*number of subdom. in y direction*/
NzClst(0), /*number of subdom. in z direction*/

nx(0), /*number of elem. in x direction*/
ny(0), /*number of elem. in y direction*/
nz(0), /*number of elem. in z direction*/

nxSub(0), /*number of elem. in x direction on sub*/
nySub(0), /*number of elem. in y direction on sub*/
nzSub(0), /*number of elem. in z direction on sub*/

cx(0), /*refinement - x*/
cy(0), /*refinement - y*/
cz(0), /*refinement - z*/

n_elementsAll(0), /*number of all elements*/
n_elementsSub(NULL),
n_nodsAll(0), /*number of all coordinates*/
n_nodsSub(NULL),
n_neighbClst(0),
n_facesAll(0), /*number of all edges*/
n_facesSub(NULL),
neqAll(0), /*number of all equations*/
neqSub(NULL),
n_subdomains(0), /*number of subdomains*/
n_clusters(0), /*number of subdomains*/
n_subdomOnClust(0), /*number of subdomains*/
n_exterDOFs(NULL),
n_exterDOFsClust(0),
fz_total(0),
flag_contact(0),
flag_redund_lagr_mult(0),
flag_linear_system(0),
flag_store_VTK(0),
flag_plotSol(0),
flag_storeCompleteMatrices(0),
flag_saveAllDataSendingToFortran(0),
vtk_min_rank(0),
vtk_max_rank(0),
verbose(0),
MPIrank(0),

eps0(0),
eps1(0),
maxIter(0),

i_load_step(0),
i_sub_step(0),
n_sub_steps(0),
n_load_steps(0),
del_deltaVec(0),
deltaVec(NULL),
vec_neqSub(NULL),
vec_neqClst(NULL),
vec_n_elementsSub(NULL),
vec_n_elementsClst(NULL),
vec_n_nodsSub(NULL),
vec_n_nodsClst(NULL),
vec_numelSub(NULL),
vec_numelClst(NULL),
vec_globalSubNumbering(NULL),
vec_globalSubNumberingMeshGen(NULL),
vec_localSubNumbering(NULL),


youngsModulus(0),
poissonsRatio(0),
density(0)
{
	//double acceleration[3];
 memset(acceleration, 0, 3 * sizeof(double));
	}

CDomain::~CDomain() {
	if (deltaVec)               delete [] deltaVec;
  if (vec_neqSub)             delete [] vec_neqSub;
  if (vec_neqClst)            delete [] vec_neqClst;
  if (vec_n_elementsSub)      delete [] vec_n_elementsSub;
  if (vec_n_elementsClst)     delete [] vec_n_elementsClst;
  if (vec_n_nodsSub)          delete [] vec_n_nodsSub;
  if (vec_n_nodsClst)         delete [] vec_n_nodsClst;
  if (vec_numelSub)           delete [] vec_numelSub;
  if (vec_numelClst)          delete [] vec_numelClst;
  if (vec_globalSubNumbering) delete [] vec_globalSubNumbering;
  if (vec_globalSubNumberingMeshGen) {
    delete [] vec_globalSubNumberingMeshGen;
  }
  if (vec_localSubNumbering)  delete [] vec_localSubNumbering;
  if (n_elementsSub)          delete [] n_elementsSub;
  if (n_nodsSub)              delete [] n_nodsSub;
  if (n_facesSub)             delete [] n_facesSub;
  if (neqSub)                 delete [] neqSub;
  if (n_exterDOFs)            delete [] n_exterDOFs;
  if (n_nodOnEdg)             delete [] n_nodOnEdg;
}


void CDomain::setFromOptions(int* argc, char*** argv) {

  int rank,MPIsize;
  MPI_Comm_rank(this->comm,&rank);
  MPI_Comm_size(comm, &MPIsize);

  int Cxyz[3] = {1,1,1}; // number of clust. in x, y, z direction					Cxyz[3] = {1,1,2};
  int Nxyz[3] = {2,2,2}; // number of subs in x, y, z direction 					Nxyz[3] = {1,2,1};
  int nxyz[3] = {10,10,10}; // number of elements in x,y,z direction (each subs.)		nxyz[3] = {1,1,1};
  
  bool nxyzset[3] = {0,0,0};
  bool Cxyzset[3] = {0,0,0};
  
  // HFETI corner setup 
  nCorners_X			= 1;
  nCorners_Y			= 1;
  nCorners_Z			= 1; 
  flag_DP_inner			= false;
  flag_DP_eges			= true;
  bool flag_contact  	= false;
  int flag_store_VTK	= 3;
  int verbose = 0;
	char c;
  int n_load_steps=1;
  int n_sub_steps=0;



  vtk_min_rank    = 0;
  vtk_max_rank    = MPIsize;
#ifndef WIN32

  nCorners_X		= 0;
  nCorners_Y		= 0;
  nCorners_Z		= 0; 
  flag_DP_inner		= false;
  flag_DP_eges		= false;
  flag_contact  	= false;
  flag_store_VTK	= false;


  int opt= 0;

  static struct option long_options[] = {
	  {"Nx",          required_argument, 0,  'a' },
	  {"Ny",          required_argument, 0,  'b' },
	  {"Nz",          required_argument, 0,  'c' },
	  {"nx",          required_argument, 0,  'd' },
	  {"ny",          required_argument, 0,  'e' },
	  {"nz",          required_argument, 0,  'f' },
	  {"Clx",         required_argument, 0,  'g' },
	  {"Cly",         required_argument, 0,  'h' },
	  {"Clz",         required_argument, 0,  'i' },
	  {"Cox",         required_argument, 0,  'j' },
	  {"Coy",         required_argument, 0,  'k' },
	  {"Coz",         required_argument, 0,  'l' },
	  {"EDGE_COR",    no_argument,       0,  'm' },
	  {"FACE_COR",    no_argument,       0,  'n' },
	  {"VTK",         required_argument, 0,  'o' },
      {"vtk_max_rank",required_argument, 0,  'p' },
      {"vtk_min_rank",required_argument, 0,  'q' },
	  {"HFETI",       required_argument, 0,  'A' },
	  {"KINV",        no_argument,       0,  'B' },
	  {"PIPECG",      required_argument, 0,  'C' },
	  {"GGTINV",      required_argument, 0,  'D' },
	  {"ITERS",       required_argument, 0,  'E' },
	  {"TOL",         required_argument, 0,  'F' },
	  {0,             0,                 0,   0  }
  };


  int long_index =0;
  while ((opt = getopt_long(*argc, *argv,"a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:A:B:C:D:E:F", 
	  long_options, &long_index )) != -1) {
		  switch (opt) {
		  case 'a' : Nxyz[0] = atoi(optarg);
			  break;
		  case 'b' : Nxyz[1] = atoi(optarg);
			  break;
		  case 'c' : Nxyz[2] = atoi(optarg);
			  break;
		  case 'd' : nxyz[0] = atoi(optarg);
			  break;
		  case 'e' : nxyz[1] = atoi(optarg); 
			  break;
		  case 'f' : nxyz[2] = atoi(optarg);
			  break;
		  case 'g' : Cxyz[0] = atoi(optarg);
			  break; 
		  case 'h' : Cxyz[1] = atoi(optarg);
			  break;
		  case 'i' : Cxyz[2] = atoi(optarg);
			  break;
		  case 'j' : nCorners_X = atoi(optarg);
			  break;
		  case 'k' : nCorners_Y = atoi(optarg);
			  break;
		  case 'l' : nCorners_Z = atoi(optarg);
			  break;
		  case 'm' : flag_DP_eges  = true;
			  break;
		  case 'n' : flag_DP_inner = true;
			  break;
		  case 'o' : flag_store_VTK = atoi(optarg);
			  break;
		  case 'p' : vtk_max_rank = atoi(optarg);
			  break;
		  case 'q' : vtk_min_rank = atoi(optarg);
			  break;
		  case 'A' : ;//cluster.USE_HFETI = atoi(optarg);
			  break;
		  case 'B' : ;//cluster.USE_KINV  = 1;
			  break;
		  case 'C' : ;//solver.USE_PIPECG = atoi(optarg);
			  break;
		  case 'D' : ;//solver.USE_GGtINV = atoi(optarg);
			  break;
		  case 'E' : ;//solver.USE_GGtINV = atoi(optarg);
			  break;
		  case 'F' : ;//solver.epsilon     = atoi(optarg);
			  break;
		  default: 
			  goto loopend2;
		  }
  }
  loopend2:


   //////  opterr = 0;
		 //////while ((c = getopt(*argc, *argv, "C:P:Q:R:x:y:z:X:Y:Z:l:s:v:cV")) != -1)
			//////	{
			////// switch (c)
			//////	 {
			//////	 case 'C':
			//////		Nxyz[0] = atoi(optarg);
			//////		Nxyz[1] = Nxyz[0];
			//////		Nxyz[2] = Nxyz[0]; 
			//////		break; 
			//////	 case 'P':
			//////		Nxyz[0] = atoi(optarg);
			//////		break; 
			//////	 case 'Q':
			//////		Nxyz[1] = atoi(optarg);
			//////		break; 
			//////	 case 'R':
			//////		Nxyz[2] = atoi(optarg);
			//////		break; 
			//////     case 'x':
			//////		 nxyz[0] = atoi(optarg);
   //////                  nxyzset[0] = 1;
			//////		 break;
			//////	 case 'y':
			//////		 nxyz[1] = atoi(optarg);
			//////		 nxyzset[1] = 1;
			//////		 break;
			//////	 case 'z':
			//////		 nxyz[2] = atoi(optarg);
			//////		 nxyzset[2] = 1;
			//////		 break;
			//////	 case 'X':
			//////		 Cxyz[0] = atoi(optarg);
			//////		 Cxyzset[0] = 1;
			//////		 break;
			//////	 case 'Y':
			//////		 Cxyz[1] = atoi(optarg);
			//////		 Cxyzset[1] = 1;
			//////		 break;
			//////	 case 'Z':
			//////		 Cxyz[2] = atoi(optarg);
			//////		 Cxyzset[2] = 1;
			//////		 break;
			//////	 case 'v':
			//////		 verbose = atoi(optarg);
			//////		 break;
			//////	 case 'l':
			//////		 n_load_steps = atoi(optarg);
			//////		 break;
			//////	 case 's':
			//////		 n_sub_steps = atoi(optarg);
			//////		 break;
			//////	 case 'c':
			//////		 flag_contact = true;
			//////		 break;
			//////	 case 'V':
			//////		 flag_store_VTK = true;//atoi(optarg)>0;
			//////		 break;
			//////	default:
   //////        goto loopend;
			//////	 }
		 ////// }

  //loopend:
  if (nxyzset[0] && !nxyzset[1]) nxyz[1] = nxyz[0];
  if (nxyzset[0] && !nxyzset[2]) nxyz[2] = nxyz[0];
  if (Cxyzset[0] && !Cxyzset[1]) Cxyz[1] = Cxyz[0];
  if (Cxyzset[0] && !Cxyzset[2]) Cxyz[2] = Cxyz[0];

#endif

	//  OPTIONS - solver, data saving,plotingSolution

//  // face 001
//  this->neuCondOnFaces.push_back(CBoundaryCondOnFace());
//  this->neuCondOnFaces[0].iFace = 1;
//  this->neuCondOnFaces[0].val_xyz[0]=0.;
//  this->neuCondOnFaces[0].val_xyz[1]=0.;
//  this->neuCondOnFaces[0].val_xyz[2]=0.;
//  // face 002
//  this->neuCondOnFaces.push_back(CBoundaryCondOnFace());
//  this->neuCondOnFaces[1].iFace = 2;
//  this->neuCondOnFaces[1].val_xyz[0]=0.;
//  this->neuCondOnFaces[1].val_xyz[1]=0.;
//  this->neuCondOnFaces[1].val_xyz[2]=0.;

  bool flag_redund_lagr_mult = true;
	bool flag_storeCompleteMatrices = false;
	bool flag_plotSol = false;
	bool flag_saveAllDataSendingToFortran = false;
	//  own solver setting  
  int  max_nDOFs_u_is_stored = 2e9; //1e5;
	double epsCG = 1.0e-7;
	double epsNewRaph = 1.0e-5;
	int maxIter = 1000;
	//
	double cxyz[3] = {1,1,1};
	double Lxyz[3] = {1., 1., 1.}; // length of cube edges
	double fz_total = -1000.; // force loaded upper face
	double acceleration[3] = {0.,0.,-0*9810.};
	double density = 7.850e-3;
	double youngsModulus = 2.0e5;
	double poissonsRatio = 0.33;

//
  int nxyzSub[3];
  int NxyzClst[3];
//		
	for (int i = 0; i < 3; i++){
    nxyzSub[i] = nxyz[i];
    NxyzClst[i] = Nxyz[i];
    Nxyz[i] = Nxyz[i]*Cxyz[i];
		nxyz[i] = nxyz[i]*Nxyz[i];
  }
	int n_subdomains    = Nxyz[0]*Nxyz[1]*Nxyz[2];
	int n_clusters      = Cxyz[0]*Cxyz[1]*Cxyz[2];
	int n_elementsAll = nxyz[0]*nxyz[1]*nxyz[2];
	int n_nodsAll = (nxyz[0] + 1)*(nxyz[1] + 1)*(nxyz[2] + 1);
	int n_facesAll =
		2 * (nxyz[0] * nxyz[1] +  nxyz[0] * nxyz[2] + nxyz[1] * nxyz[2]);
//-+-
  int n_subdomOnClust = NxyzClst[0]*NxyzClst[1]*NxyzClst[2];
  int *n_elementsSub = new int[n_subdomOnClust];
  int *n_nodsSub     = new int[n_subdomOnClust];
  int *n_facesSub    = new int[n_subdomOnClust]; 
  for (int i = 0;i<n_subdomOnClust;i++){
    n_elementsSub[i] = nxyzSub[0] * nxyzSub[1] * nxyzSub[2];
    n_nodsSub[i]     = (nxyzSub[0]+1)*(nxyzSub[1]+1)*(nxyzSub[2]+1);
    n_facesSub[i]    = 2 * (nxyzSub[0] * nxyzSub[1] +  
                            nxyzSub[0] * nxyzSub[2] + 
                            nxyzSub[1] * nxyzSub[2]);
  }
//
  int *neqSub       = new int[n_subdomOnClust];
  int *n_exterDOFs  = new int[n_subdomOnClust];
  int *n_nodOnEdg   = new int[n_subdomOnClust];
  this->vec_globalSubNumbering = new int[n_subdomOnClust];
  this->vec_globalSubNumberingMeshGen = new int[n_subdomOnClust];
  this->vec_localSubNumbering = new int[n_subdomOnClust];
//-+-
 // ................................................................. 

  for (int i=0;i<n_subdomOnClust;i++){
    this->vec_globalSubNumbering[i] = rank*n_subdomOnClust + i;
    this->vec_localSubNumbering[i] = i;
  }




	this->Lx      = Lxyz[0];
	this->Ly      = Lxyz[1];
	this->Lz      = Lxyz[2];
	this->Nx      = Nxyz[0];
	this->Ny      = Nxyz[1];
	this->Nz      = Nxyz[2];
	this->Cx      = Cxyz[0];
	this->Cy      = Cxyz[1];
	this->Cz      = Cxyz[2];
	this->NxClst  = NxyzClst[0];
	this->NyClst  = NxyzClst[1];
	this->NzClst  = NxyzClst[2];
	this->nx      = nxyz[0];
	this->ny      = nxyz[1];
	this->nz      = nxyz[2];
	this->nxSub   = nxyzSub[0];
	this->nySub   = nxyzSub[1];
	this->nzSub   = nxyzSub[2];
	this->cx      = cxyz[0];
	this->cy      = cxyz[1];
	this->cz      = cxyz[2];
	this->n_nodsAll = n_nodsAll;
  this->n_nodsSub = n_nodsSub;
	this->n_elementsSub= n_elementsSub;
	this->n_elementsAll = n_elementsAll;
	this->n_subdomOnClust = n_subdomOnClust;
	this->n_clusters= n_clusters;
	this->n_facesAll = n_facesAll;
	this->n_facesSub = n_facesSub;
	this->neqAll = n_nodsAll * 3;
  this->deltaVec = new double[n_load_steps];
	this->n_subdomains = n_subdomains;
	//
  this->flag_contact = flag_contact;
  gv_flag_linear_system = (n_load_steps==1 && n_sub_steps==0);
  //this->flag_linear_system = gv_flag_linear_system;
	this->flag_store_VTK = flag_store_VTK;
	this->flag_plotSol = flag_plotSol;
	this->flag_storeCompleteMatrices = flag_storeCompleteMatrices;
	this->flag_redund_lagr_mult= flag_redund_lagr_mult;
	this->verbose = verbose;
  this->MPIrank = rank;
  this->max_nDOFs_u_is_stored = max_nDOFs_u_is_stored;
	//
	this->eps0 = epsCG;
	this->eps1 = epsNewRaph;
	this->maxIter = maxIter;
	this->n_load_steps = n_load_steps;
	this->n_sub_steps = n_sub_steps;
	this->fz_total = fz_total;
	this->density = density;
	this->youngsModulus = youngsModulus;
	this->poissonsRatio = poissonsRatio;
	this->flag_saveAllDataSendingToFortran = flag_saveAllDataSendingToFortran;
  this->neqSub       = neqSub;
  this->n_exterDOFs  = n_exterDOFs;
  this->n_nodOnEdg   = n_nodOnEdg;
	//
	for (int i = 0;i<3;i++){
		this->acceleration[i] = acceleration[i];
	}
}

void CDomain::print(){
  int rank;
  MPI_Comm_rank(this->comm,&rank);
  if (rank) return;
  int nxs=this->nxSub;
  int nys=this->nySub;
  int nzs=this->nzSub;
  
  printf("----------------------------------------------------------------------\n");
  printf("#subdomains         %16d  (Nx Ny Nz)=(%d %d %d)\n", 
                          this->n_subdomains, this->Nx, this->Ny, this->Nz);
  printf("#subdomains on cls. %16d             (%d %d %d)\n", 
                          this->n_subdomOnClust, this->NxClst, this->NyClst, this->NzClst);
  printf("#clusters           %16d  (Cx Cy Cz)=(%d %d %d)\n",  
                          this->n_clusters, this->Cx, this->Cy, this->Cz);
  printf("#elements per subd. %16d  (nx ny nz)=(%d %d %d)\n", 
                          this->n_elementsSub[0],nxs,nys,nzs);
  printf("#all elements:      %16d  (nx ny nz)=(%d %d %d)\n",   
                          this->n_elementsAll, this->nx, this->ny, this->nz); printf("#coordinates:       %16d\n", this->n_nodsAll);
  printf("#DOFs undecomp.:    %16d\n", this->neqAll);
  printf("#DOFs decomp.:      %16d\n", 3*(nxs+1)*(nys+1)*(nzs+1)*this->n_subdomains);
  printf("contact:            %16d\n", this->flag_contact);
  printf("----------------------------------------------------------------------\n");
}
