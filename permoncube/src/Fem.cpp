#include "Fem.h"


CFem::CFem(MPI_Comm comm, int _i_domOnClust) {
  this->comm = comm;
  i_domOnClust = _i_domOnClust;
}

CFem::CFem(const CFem& other):
  comm(other.comm),
  mesh(other.mesh),
  bound_cond(other.bound_cond),
  i_domOnClust(other.i_domOnClust)
{

}

CFem::~CFem() {
}

void CFem::coordinate_generator3d(Coordinates &coordinates)
{
	int clusterCount[3] = { 1, 1, 1 };
	int subInClusterCount[3] = { 2, 2, 2 };
	int elementsInSubdomainCount[3] = { 10, 10, 10 };

	int nnx = subInClusterCount[0] * elementsInSubdomainCount[0] + 1;
	int nny = subInClusterCount[1] * elementsInSubdomainCount[1] + 1;
	int nnz = subInClusterCount[2] * elementsInSubdomainCount[2] + 1;
	double lenght[3] = { 1, 1, 1 };

	double stepx = lenght[0] / (nnx - 1);
	double stepy = lenght[1] / (nny - 1);
	double stepz = lenght[2] / (nnz - 1);
	idx_t index = 0;
	for (int z = 0; z < nnz; z++) {
		for (int y = 0; y < nny; y++) {
			for (int x = 0; x < nnx; x++) {
				coordinates[index++] = Point(x * stepx, y * stepy, z * stepz);
			}
		}
	}
}

void CFem::element_generator3d(Mesh &mesh)
{

	//	###################################################
	//	#                                                 #
	//	#             A z-coord.                          #
	//	#             |                                   #
	//	#             |            E3                     #
	//	#             |_ _ _ _ _ _ _                      #
	//	#            /     E5      /|                     #
	//	#           /_ _ _ _ _ _  / |                     #
	//	#          |      |      |  |                     #
	//	#        E4|      |      |E2|                     #
	//	#          |_ _ _ |_ _ _ |  |       y-coord.      #
	//	#          |    E1|      |  |------->             #
	//	#          |      |      | /                      #
	//	#          |_ _ _ |_ _ _ |/                       #
	//	#         /                                       #
	//	#        /       E0                               #
	//	#       /                                         #
	//	#      v  x-coord.                                #
	//	#                                                 #
	//	###################################################


	int clusterCount[3] = { 1, 1, 1 };
	int subInClusterCount[3] = { 2, 2, 2 };
	int elementsInSubdomainCount[3] = { 10, 10, 10 };

	idx_t indices[8];
	int nnx = subInClusterCount[0] * elementsInSubdomainCount[0] + 1;
	int nny = subInClusterCount[1] * elementsInSubdomainCount[1] + 1;
	int nnz = subInClusterCount[2] * elementsInSubdomainCount[2] + 1;

	int offset[3];

	for (int subz = 0; subz < subInClusterCount[2]; subz++) {
		for (int suby = 0; suby < subInClusterCount[1]; suby++) {
			for (int subx = 0; subx < subInClusterCount[0]; subx++) {
				offset[2] = subz * elementsInSubdomainCount[2];
				offset[1] = suby * elementsInSubdomainCount[1];
				offset[0] = subx * elementsInSubdomainCount[0];
				for (int z = offset[2]; z < offset[2] + elementsInSubdomainCount[2]; z++) {
					for (int y = offset[1]; y < offset[1] + elementsInSubdomainCount[1]; y++) {
						for (int x = offset[0]; x < offset[0] + elementsInSubdomainCount[0]; x++) {
							indices[0] = nnx * nny *  z      + nnx *  y      + x;
							indices[1] = nnx * nny *  z      + nnx *  y      + x + 1;
							indices[2] = nnx * nny *  z      + nnx * (y + 1) + x + 1;
							indices[3] = nnx * nny *  z      + nnx * (y + 1) + x;
							indices[4] = nnx * nny * (z + 1) + nnx *  y      + x;
							indices[5] = nnx * nny * (z + 1) + nnx *  y      + x + 1;
							indices[6] = nnx * nny * (z + 1) + nnx * (y + 1) + x + 1;
							indices[7] = nnx * nny * (z + 1) + nnx * (y + 1) + x;
							mesh.push_element(new Hexahedron(indices));
						}
					}
				}

			}
		}
	}
}

void CFem::mesh_generator3d(CDomain *domainG) {



//	###################################################
//	#                                                 #
//	#             A z-coord.                          #
//	#             |                                   #
//	#             |            E3                     #
//	#             |_ _ _ _ _ _ _                      #
//	#            /     E5      /|                     #
//	#           /_ _ _ _ _ _  / |                     #
//	#          |      |      |  |                     #
//	#        E4|      |      |E2|                     #
//	#          |_ _ _ |_ _ _ |  |       y-coord.      #
//	#          |    E1|      |  |------->             #
//	#          |      |      | /                      #
//	#          |_ _ _ |_ _ _ |/                       #
//	#         /                                       #
//	#        /       E0                               #
//	#       /                                         #
//	#      v  x-coord.                                #
//	#                                                 #
//	###################################################

	int	globoalIndSub = this->mesh.i_domGlobalMeshGen;


	int nx = domainG->nx;
	int ny = domainG->ny;
	int nz = domainG->nz;

	int Nx = domainG->Nx;
	int Ny = domainG->Ny;
	int Nz = domainG->Nz;

	double cx = domainG->cx;
	double cy = domainG->cy;
	double cz = domainG->cz;

	double hx = 1.0 / nx;
	double hy = 1.0 / ny;
	double hz = 1.0 / nz;

	int nnx = nx + 1;
	int nny = ny + 1;
	int nnxy = nnx * nny;
	int cnt = 0;
	int m = 0;
	CSolid45 * element = this->mesh.element;
	CRectangularFace * edgeSub = this->mesh.faceSub;


  int nxSub = domainG->nxSub; 
  int nySub = domainG->nySub; 
  int nzSub = domainG->nzSub;
  int NxClst= domainG->NxClst; 
  int NyClst= domainG->NyClst; 
  int NzClst= domainG->NzClst;



  // + + + only for the cube - start
  int Isub,Jsub,Ksub;
  Ksub = ceil(double (globoalIndSub+1)/(Nx*Ny));
  int inXYplane = (globoalIndSub+1)-(Ksub-1)*Nx*Ny;
  Jsub = ceil(double( inXYplane)/Nx);
  Isub = inXYplane - Nx*(Jsub-1);
  // + + + only for the cube - end

	cnt = 0;

  m = 0;
  Isub = Isub-1;
  Jsub = Jsub-1;
  Ksub = Ksub-1;
  int m1=0;
	for (int k = 0; k < nzSub; k++) {
		for (int j = 0; j < nySub; j++) {
			for (int i = 0; i < nxSub; i++) {
				element[m].inod_glob[0] = 0 + 
							(i+Isub*nxSub) + (j+Jsub*nySub) * nnx + (k+Ksub*nzSub) * nnxy;
				element[m].inod_glob[1] = 1 + 
							(i+Isub*nxSub) + (j+Jsub*nySub) * nnx + (k+Ksub*nzSub) * nnxy;
				element[m].inod_glob[2] = 1 + 
							(i+Isub*nxSub) + (j+Jsub*nySub + 1) * nnx + (k+Ksub*nzSub) * nnxy;
				element[m].inod_glob[3] = 0 + 
							(i+Isub*nxSub) + (j+Jsub*nySub + 1) * nnx + (k+Ksub*nzSub) * nnxy;
				element[m].inod_glob[4] = 0 + 
							(i+Isub*nxSub) + (j+Jsub*nySub) * nnx + (k+Ksub*nzSub + 1) * nnxy;
				element[m].inod_glob[5] = 1 + 
							(i+Isub*nxSub) + (j+Jsub*nySub) * nnx + (k+Ksub*nzSub + 1) * nnxy;
				element[m].inod_glob[6] = 1 + 
							(i+Isub*nxSub) + (j+Jsub*nySub + 1) * nnx + (k+Ksub*nzSub + 1) * nnxy;
				element[m].inod_glob[7] = 0 + 
							(i+Isub*nxSub) + (j+Jsub*nySub + 1) * nnx + (k+Ksub*nzSub + 1) * nnxy;
				element[m].ordinalNumber = m;

				// face 0
        if (k == 0) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[0];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[3];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[2];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[1];
					edgeSub[m1].iFaceSub = 0;
  				edgeSub[m1].iFace= -1;
					if ((k+Ksub*nzSub)==0){
  					edgeSub[m1].iFace= 0;
					}
					edgeSub[m1].iElem = m;
  				m1++;
				}
				// face 1
        if (i == nxSub-1) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[1];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[2];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[6];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[5];
					edgeSub[m1].iFaceSub = 1;
  				edgeSub[m1].iFace= -1;
					if ((i+Isub*nxSub)==(nx-1)){
  					edgeSub[m1].iFace= 1;
					}
					edgeSub[m1].iElem = m;
				  m1++;
				}
				// face 2
        if (j == nySub-1) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[7];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[6];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[2];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[3];
					edgeSub[m1].iFaceSub = 2;
  				edgeSub[m1].iFace= -1;
					if ((j + Jsub*nySub) == (ny-1)){
  					edgeSub[m1].iFace= 2;
					}
					edgeSub[m1].iElem = m;
				  m1++;
				}
				// face 3
        if (i == 0) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[4];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[7];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[3];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[0];
					edgeSub[m1].iFaceSub = 3;
  				edgeSub[m1].iFace= -1;
					if ((i+Isub*nxSub)==0){
  					edgeSub[m1].iFace= 3;
					}
					edgeSub[m1].iElem = m;
				  m1++;
				}
				// face 4
        if (j == 0) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[0];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[1];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[5];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[4];
					edgeSub[m1].iFaceSub = 4;
  					edgeSub[m1].iFace= -1;
					if ((j+Jsub*nySub)==0){
  					edgeSub[m1].iFace= 4;
					}
					edgeSub[m1].iElem = m;
  				m1++;
				}
				// face 5
        if (k == nzSub-1) {
					edgeSub[m1].inod_glob[0] = element[m].inod_glob[4];
					edgeSub[m1].inod_glob[1] = element[m].inod_glob[5];
					edgeSub[m1].inod_glob[2] = element[m].inod_glob[6];
					edgeSub[m1].inod_glob[3] = element[m].inod_glob[7];
					edgeSub[m1].iFaceSub = 5;
  				edgeSub[m1].iFace= -1;
					if ((k+Ksub*nzSub)==(nz-1)){
  					edgeSub[m1].iFace= 5;
					}
					edgeSub[m1].iElem = m;
  				m1++;
				}
				m++;
			}//i
		}//j
	}//k



  cnt=0;

  Isub = Isub+1;
  Jsub = Jsub+1;
  Ksub = Ksub+1;




//-----------------------------------------------------------------
  std::map<longInt,int> & g2l_nodes = this->mesh.g2l_nodes;
  std::vector<longInt> & l2g_nodes = this->mesh.l2g_nodes;

  int cntijv=0;
  int n_nodFromAllSubEl= domainG->n_elementsSub[this->i_domOnClust]*8;
  l2g_nodes.resize(n_nodFromAllSubEl);

  cnt=0;

  for (int i = 0; i < domainG->n_elementsSub[this->i_domOnClust]; i++) {
		for (int j = 0; j < 8; j++) {
			l2g_nodes[cntijv] = element[cnt].inod_glob[j];
			cntijv++;
		} // loop j
    cnt++;
  } // loop i
  //
  qsort(&(l2g_nodes[0]),n_nodFromAllSubEl, sizeof (longInt), CLinearAlgebra::compareLongInt);
  std::vector<longInt>::iterator it;
  it = std::unique (l2g_nodes.begin(), l2g_nodes.end());
  l2g_nodes.resize( std::distance(l2g_nodes.begin(),it) );
  for (unsigned int ii = 0; ii < l2g_nodes.size(); ii++){
    g2l_nodes.insert ( std::pair<longInt,int>(l2g_nodes[ii],ii) );
  }
  //
  // optimized version of generation of coordinates over one subdomain

  Ksub--; Jsub--; Isub--;
  int m0 = 0;
	for (int k = Ksub*nzSub; k < (Ksub+1)*nzSub + 1; k++) {
	  for (int j = Jsub*nySub; j < (Jsub+1)*nySub + 1; j++) {
	    for (int i = Isub*nxSub; i < (Isub+1)*nxSub + 1; i++) {
					this->mesh.coordinateSub[m0].x = (pow(i*hx,cx))*domainG->Lx;
					this->mesh.coordinateSub[m0].y = (pow(j*hy,cy))*domainG->Ly;
					this->mesh.coordinateSub[m0].z = (pow(k*hz,cz))*domainG->Lz;
				  m0++;
			}
		}
	}
  Ksub++; Jsub++; Isub++;

//	m = 0;
//  int m0 = 0;
//    // TODO 20150224
//	// coordinates0 matrix
//	for (int k = 0; k < domainG->nz + 1; k++) {
//		for (int j = 0; j < domainG->ny + 1; j++) {
//			for (int i = 0; i < domainG->nx + 1; i++) {
//        if (m0<domainG->n_nodsSub[this->i_domOnClust] && this->mesh.l2g_nodes[m0]==m){
//					this->mesh.coordinateSub[m0].x = (pow(i*hx,cx))*domainG->Lx;
//					this->mesh.coordinateSub[m0].y = (pow(j*hy,cy))*domainG->Ly;
//					this->mesh.coordinateSub[m0].z = (pow(k*hz,cz))*domainG->Lz;
//				  m0++;
//				}
//				m++;
//			}
//		}
//	}
  
// Finding of external nodes on cube defined by (nx,ny,nz) elements
//					(the same process on all ranks)
  int nnxSub = nxSub+1;
  int nnySub = nySub+1;
  int nnzSub = nzSub+1;

  int cntN = 0;
  cnt = 0;
  for (int k = 0;k<nnzSub;k++){
		for (int j = 0;j<nnySub;j++){
			for (int i = 0;i<nnxSub;i++){
        if (k==0||k==nnzSub-1){
				  this->mesh.nodOnEdg[cntN]=cnt; 
				  cntN++;
				}
				else {
          if ((i==0 || i==nnxSub-1 || j==0 || j==nnySub-1)){
				  	this->mesh.nodOnEdg[cntN]=cnt; 
						cntN++;
					}
				}
				cnt++;
			}
		}
	}


// A2,start: Finding neigbours subdomains of each subdomain
  int *neighbSub = new int[26];
  cnt = 0;
  int tmpRankGuessed;
//

  for (int k = Ksub-1;k<=Ksub+1;k++){
		for (int j = Jsub-1;j<=Jsub+1;j++){
		  for (int i = Isub-1;i<=Isub+1;i++){
				tmpRankGuessed = Nx*Ny*(k-1)+Nx*(j-1)+i-1;
        if ((i>0 && j>0 && k>0) && 
						(i<Nx+1 && j<Ny+1 && k<Nz+1) &&
						(tmpRankGuessed != globoalIndSub)) {
			     neighbSub[cnt] =	tmpRankGuessed;
           cnt++;
				}
		  }
		}
  }

  domainG->n_neighbSub = cnt;
  this->mesh.neighbSub = neighbSub;
  int nnrow = (nxSub+1)*(nySub+1)*nzSub;
 // 
  int *fixingNodes= new int[8];
  fixingNodes[0] = 0;
  fixingNodes[1] = nxSub;
  fixingNodes[2] = (nxSub+1)*(nySub+1)-1-nxSub;
  fixingNodes[3] = (nxSub+1)*(nySub+1)-1;    
  fixingNodes[4] = 0 + nnrow;
  fixingNodes[5] = nxSub + nnrow;
  fixingNodes[6] = (nxSub+1)*(nySub+1)-1-nxSub + nnrow;
  fixingNodes[7] = (nxSub+1)*(nySub+1)-1 + (nxSub+1)*(nySub+1)*nzSub;    
 // 
  int *fixingDOFs = new int[24];
 // 
  for (int i = 0;i<8;i++){
    fixingDOFs[3*i+0] = 3*fixingNodes[i]+0;
    fixingDOFs[3*i+1] = 3*fixingNodes[i]+1;
    fixingDOFs[3*i+2] = 3*fixingNodes[i]+2;
  }
 // 
  this->mesh.fixingDOFs = fixingDOFs;
  this->mesh.fixingNodes= fixingNodes;
//


//  vector <int> myVec;
  int n_myVec;
  n_myVec=(nxSub+1)*(nySub+1)*(nzSub+1)-(nxSub-1)*(nySub-1)*(nzSub-1);
//  myVec.reserve(n_myVec);
  this->mesh.DP_NodesAll.reserve(n_myVec);
  vector<int> tmpVec;
  tmpVec.resize(8);

 // number of 'corners' per edge 

  longInt tmp_i;
  longInt I,J,K;
  longInt remI,remJ,remK;
  longInt tmpArray[8];
  int cnt_tmpArray;



//TODO 1) nCorners will be specified by user (extension of getopt)
//TODO 2) collection  * flag_DP_inner, 
//                    * flag_DP_eges and  
//                    * flag_DP_real_corners) 
//        add to getopt
//  
  
  bool flag_DP_inner = domainG->flag_DP_inner; 
  bool flag_DP_eges  = domainG->flag_DP_eges;
  int nCorners_X     = domainG->nCorners_X; // 3;
  int nCorners_Y     = domainG->nCorners_Y; // 3;
  int nCorners_Z     = domainG->nCorners_Z; // 3;
// correction leading to considered number of 'corners' points per edge(s)
  nCorners_X+=1;
  nCorners_Y+=1;
  nCorners_Z+=1;

  int stepX=max(nxSub/nCorners_X,1);
  int stepY=max(nySub/nCorners_Y,1);
  int stepZ=max(nzSub/nCorners_Z,1);
  bool boo_cond1,boo_cond2,boo_cond1_in,boo_cond2_in;
//
  for (int i=0;i<=int(nxSub/2);i+=stepX){
    for (int j=0;j<=int(nySub/2);j+=stepY){
      for (int k=0;k<=int(nzSub/2);k+=stepZ){
//

        boo_cond1=flag_DP_inner && (i==0 || j==0 || k==0);
        boo_cond2=flag_DP_eges && ((int(i==0)+int(j==0)+int(k==0))>1);
//
        if (boo_cond1 || boo_cond2){
          tmpVec.resize(8);
          tmpVec[0] =      i +(nxSub+1)*     j  +      k  * (nxSub+1)*(nySub+1);
          tmpVec[1] = (nxSub-i)+(nxSub+1)*     j  +      k  * (nxSub+1)*(nySub+1);
          tmpVec[2] =      i +(nxSub+1)*(nySub-j) +      k  * (nxSub+1)*(nySub+1);
          tmpVec[3] = (nxSub-i)+(nxSub+1)*(nySub-j) +      k  * (nxSub+1)*(nySub+1);
          tmpVec[4] =      i +(nxSub+1)*     j  + (nzSub-k) * (nxSub+1)*(nySub+1);
          tmpVec[5] = (nxSub-i)+(nxSub+1)*     j  + (nzSub-k) * (nxSub+1)*(nySub+1);
          tmpVec[6] =      i +(nxSub+1)*(nySub-j) + (nzSub-k) * (nxSub+1)*(nySub+1);
          tmpVec[7] = (nxSub-i)+(nxSub+1)*(nySub-j) + (nzSub-k) * (nxSub+1)*(nySub+1);
//
          qsort(&(tmpVec[0]),8, sizeof (int), CLinearAlgebra::compareInt);
          std::vector<int>::iterator it_int; 
		  it_int = unique (tmpVec.begin(), tmpVec.end());
          tmpVec.resize( distance(tmpVec.begin(),it_int));
          cnt_tmpArray=0;
          for (it_int=tmpVec.begin();it_int!=tmpVec.end();it_int++){
            tmp_i = this->mesh.l2g_nodes[*it_int];
            K = ceil(double (tmp_i+1)/((nx+1)*(ny+1)));
            int inXYplane = (tmp_i+1)-(K-1)*(nx+1)*(ny+1);
            J = ceil(double( inXYplane)/(nx+1));
            I = inXYplane - (nx+1)*(J-1);
//
            I--; J--; K--;
            remI=I%(nxSub*NxClst);
            remJ=J%(nySub*NyClst);
            remK=K%(nzSub*NzClst);
//

            boo_cond1_in=flag_DP_inner && (remI!=0 && remJ!=0 && remK!=0);
            boo_cond2_in=flag_DP_eges && ((int(remI!=0)+int(remJ!=0)+int(remK!=0))>1);
//
            if ((boo_cond1 && boo_cond1_in) || (boo_cond2 && boo_cond2_in)){
              tmpArray[cnt_tmpArray]=tmp_i;
              cnt_tmpArray++;
            }
          }
          //myVec.insert(myVec.end(),tmpArray,tmpArray+cnt_tmpArray);
          this->mesh.DP_NodesAll.insert(this->mesh.DP_NodesAll.end(),tmpArray,tmpArray+cnt_tmpArray);

        }
      }
    }
  }
//
  if (this->mesh.DP_NodesAll.size() > 0) {
	  qsort(&(this->mesh.DP_NodesAll[0]),this->mesh.DP_NodesAll.size(), 
			  sizeof (longInt), CLinearAlgebra::compareLongInt);
	  it = unique (this->mesh.DP_NodesAll.begin(), this->mesh.DP_NodesAll.end());
	  this->mesh.DP_NodesAll.resize( distance(this->mesh.DP_NodesAll.begin(),it));
  }
//
//  if (this->i_domOnClust==0){
  //for ( it = this->mesh.DP_NodesAll.begin();it!= this->mesh.DP_NodesAll.end();it++){
  //    printf("dom: %d, point: %d\n",this->i_domOnClust,*it);
  //}


  //if (!globoalIndSub) {
	 // printf("mesh: \tdone\n");
  //}


//#ifndef WIN32
//  if (!globoalIndSub) { cout << " Mesh done                                                                "; system("date +%T.%6N"); }
//#endif
  }

//}
//
//
//
void CFem::dataDirBCSub(int *i_face,int n_facesWithDirBC, 
				int n_facesSub, int *dir_xyz) {
  int cnt=0;


  for (int i = 0;i<n_facesSub;i++){
    for (int j=0;j<n_facesWithDirBC;j++){
			if (this->mesh.faceSub[i].iFace==i_face[j]){
				cnt++;
			}
    }
  }  

  int nVec;
  if (cnt>0){
    std::vector<longInt> tmpVec;
		nVec = cnt*4;
		tmpVec.resize(nVec);

		cnt = 0; 
		for (int i = 0;i<n_facesSub;i++){
			for (int j=0;j<n_facesWithDirBC;j++){
				if (this->mesh.faceSub[i].iFace==i_face[j]){
					tmpVec[cnt  ] = this->mesh.faceSub[i].inod_glob[0];	
					tmpVec[cnt+1] = this->mesh.faceSub[i].inod_glob[1];	
					tmpVec[cnt+2] = this->mesh.faceSub[i].inod_glob[2];	
					tmpVec[cnt+3] = this->mesh.faceSub[i].inod_glob[3];	
					cnt+=4;
				}
			}
		}  

		qsort(&(tmpVec[0]),nVec, sizeof (longInt), CLinearAlgebra::compareLongInt);
		std::vector<longInt>::iterator it;
		it = std::unique (tmpVec.begin(), tmpVec.end());
		tmpVec.resize( std::distance(tmpVec.begin(),it));
    int n_tmpVec = tmpVec.size();
    for (int i = 0;i<n_tmpVec;i++){
       tmpVec[i] = this->mesh.g2l_nodes[tmpVec[i]];
    }
    int n_uniqVec = tmpVec.size();
		this->bound_cond.dirBCSub->n =   3*n_uniqVec;
		this->bound_cond.dirBCSub->ind = new int[3*n_uniqVec];
		this->bound_cond.dirBCSub->val = new double[3*n_uniqVec];

		for (int i = 0; i<n_uniqVec;i++) {
			this->bound_cond.dirBCSub->ind[3*i]   = (int)3*tmpVec[i];
			this->bound_cond.dirBCSub->val[3*i]   = 0.0;
			this->bound_cond.dirBCSub->ind[3*i+1] = (int)3*tmpVec[i]+1;
			this->bound_cond.dirBCSub->val[3*i+1] = 0.0;
			this->bound_cond.dirBCSub->ind[3*i+2] = (int)3*tmpVec[i]+2;
      double zShift = this->mesh.coordinateSub[tmpVec[i]].x + 
                      this->mesh.coordinateSub[tmpVec[i]].y;
			this->bound_cond.dirBCSub->val[3*i+2] = 0.0000*zShift;
		}
//#define test001
#ifdef test001
    this->bound_cond.dirBCSub->ind[0] = 0;
    this->bound_cond.dirBCSub->val[0] = 0.0;
    
    this->bound_cond.dirBCSub->ind[1] = 1;
    this->bound_cond.dirBCSub->val[1] = 0.0;

    this->bound_cond.dirBCSub->ind[2] = 2;
    this->bound_cond.dirBCSub->val[2] = 0.0;

    this->bound_cond.dirBCSub->ind[3] = 4;
    this->bound_cond.dirBCSub->val[3] = 0.0;

    this->bound_cond.dirBCSub->ind[4] = 5;
    this->bound_cond.dirBCSub->val[4] = 0.0;

    this->bound_cond.dirBCSub->ind[5] = 6;
    this->bound_cond.dirBCSub->val[5] = 0.0;

    this->bound_cond.dirBCSub->ind[6] = 8;
    this->bound_cond.dirBCSub->val[6] = 0.0;

    this->bound_cond.dirBCSub->ind[7] = 11;
    this->bound_cond.dirBCSub->val[7] = 0.0;

    this->bound_cond.dirBCSub->ind[8] = 12;
    this->bound_cond.dirBCSub->val[8] = 0.0;

    this->bound_cond.dirBCSub->ind[9] = 13;
    this->bound_cond.dirBCSub->val[9] = 0.0;

    this->bound_cond.dirBCSub->ind[10] = 16;
    this->bound_cond.dirBCSub->val[10] = 0.0;

    this->bound_cond.dirBCSub->ind[11] = 18;
    this->bound_cond.dirBCSub->val[11] = 0.0;
#endif

  }
  else 
  {
		this->bound_cond.dirBCSub->n =   0;
		this->bound_cond.dirBCSub->ind = NULL;
		this->bound_cond.dirBCSub->val = NULL;
  }

}


void CFem::dataConBCSub(CDomain *domainG,int *i_face,int n_facesWithConBC, 
				int n_facesSub, int *dir_xyz) {
//
  int cnt=0;
  int rank;
  MPI_Comm_rank(this->comm,&rank);
//
  if (domainG->flag_contact){
		for (int i = 0;i<n_facesSub;i++){
			for (int j=0;j<n_facesWithConBC;j++){
				if (this->mesh.faceSub[i].iFace==i_face[j]){
					cnt++;
				}
			}
		}  
  }

  if (cnt>0){
		std::vector<longInt> tmpVec;
		int nVec;
		nVec = cnt*4;
		tmpVec.resize(nVec);

		cnt = 0; 
		for (int i = 0;i<n_facesSub;i++){
			for (int j=0;j<n_facesWithConBC;j++){
				if (this->mesh.faceSub[i].iFace==i_face[j]){
					tmpVec[cnt  ] = this->mesh.faceSub[i].inod_glob[0];	
					tmpVec[cnt+1] = this->mesh.faceSub[i].inod_glob[1];	
					tmpVec[cnt+2] = this->mesh.faceSub[i].inod_glob[2];	
					tmpVec[cnt+3] = this->mesh.faceSub[i].inod_glob[3];	
					cnt+=4;
				}
			}
		}  

		qsort(&(tmpVec[0]),nVec, sizeof (longInt), CLinearAlgebra::compareLongInt);
		std::vector<longInt>::iterator it;
		it = std::unique (tmpVec.begin(), tmpVec.end());
		tmpVec.resize( std::distance(tmpVec.begin(),it));
    int n_tmpVec = tmpVec.size();
    for (int i = 0;i<n_tmpVec;i++){
       tmpVec[i] = this->mesh.g2l_nodes[tmpVec[i]];
    }
    int n_uniqVec = tmpVec.size();
		this->bound_cond.conBCSub->n =   n_uniqVec;
		this->bound_cond.conBCSub->ind = new int[n_uniqVec];
		this->bound_cond.conBCSub->val = new double[n_uniqVec];


		for (int i = 0; i<n_uniqVec;i++) {
			this->bound_cond.conBCSub->ind[i+0] = (int)(3*tmpVec[i]+0);
			this->bound_cond.conBCSub->val[i+0] = 1.0;
		}
    if (domainG->verbose>1) printf("[%d] PermonCube: subdomain with contact \n",rank);
  }
  else 
  {
		this->bound_cond.conBCSub->n =   0;
		this->bound_cond.conBCSub->ind = NULL;
		this->bound_cond.conBCSub->val = NULL;
  }
}


void CFem::applicationDirBCtoAxb(CKSparse *KSparse,double *f,CDomain *domainG){
  int n = domainG->neqSub[this->i_domOnClust], k=0, l=0;// MPIrank, MPIsize;
  int nDir = this->bound_cond.dirBCSub->n;
  double rho = 1.e5,AverDiagK = 0.0;
//  MPI_Comm_rank(comm, &MPIrank);
//  MPI_Comm_size(comm, &MPIsize);
  //
  for (int i = 0; i < n; i++) {
    AverDiagK += KSparse->val[KSparse->row_ptr[i]];
  }
  AverDiagK = AverDiagK / n;
  //
  for (int i=0;i<n;i++){
    l = k+0;
    if ((k<nDir) && (i==this->bound_cond.dirBCSub->ind[k])){
      for (int j=KSparse->row_ptr[i];j<KSparse->row_ptr[i+1];j++){
        if (j==KSparse->row_ptr[i]){       
          KSparse->val[j] = 1.0*AverDiagK;
          f[i] = this->bound_cond.dirBCSub->val[k]*AverDiagK;
        } // K[i,i] ... diagonal element
        else {
          if ((l<nDir) && (KSparse->col_ind[j]==this->bound_cond.dirBCSub->ind[l])){
            l++;
          }
          else{
            f[KSparse->col_ind[j]] -= 
                     KSparse->val[j]*this->bound_cond.dirBCSub->val[k];
          }
          KSparse->val[j] = 0.0;
        }
      }
      k++;
    }
    else {
      for (int j=KSparse->row_ptr[i];j<KSparse->row_ptr[i+1];j++){
        if (l<nDir && KSparse->col_ind[j]==this->bound_cond.dirBCSub->ind[l]){
          f[i] -= KSparse->val[j]*this->bound_cond.dirBCSub->val[l];
          KSparse->val[j] = 0.0;
          l++;
        }
      }
    }
  }
 
  // K, and f printed to files +++++++++++++++++++++++++++++++++++++++++++ start
  clock_t begin = clock();
  char filenameKe[128];
  sprintf(filenameKe, "data/K_%d.dat", 0);
  ofstream f_loc_stif_mat(filenameKe);
  for (int i = 0; i < domainG->neqSub[this->i_domOnClust]; i++) {
    for (int j = KSparse->row_ptr[i]; j < KSparse->row_ptr[i + 1]; j++) {
      f_loc_stif_mat << i << "\t" << KSparse->col_ind[j] << "\t"
        << setprecision(12) << KSparse->val[j] << endl;
    }
  }
  f_loc_stif_mat.close(); //                                000
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  ofstream frhs_f("data/f_0.dat");
  for (int i = 0; i < n; i++) {
    frhs_f << f[i] << endl;
  }
  frhs_f.close(); //--
  // K, and f printed to files +++++++++++++++++++++++++++++++++++++++++++ end
}



void CFem::dataNeuBC(double * f, int n, double fz_total,int i_face, int n_faces) {
//																																			  *
//	Neumann BC applicated on face  i_face																	*
//	        _____                                            						  *
//	     |>|_|_|_| --->                                      						  *
//	 Dir |>|_|_|_| --->  Neu                                 						  *
//	     |>|_|_|_| --->                                      						  *
//	                                                         						  *
//	    y ^                                                  						  *
//	      |                                                  						  *
//	      |                                                  						  *
//	    x o----> z                                           						  *
//	
  for (int i = 0;i<n;i++){
    f[i] = 0.0; 
  }
  
	int iDof = 0;
	double AreaOfEl,sArea=0.0;
	int in0,in1,in2;

	for (int i = 0;i<n_faces;i++){
		if (this->mesh.faceSub[i].iFace==i_face){
			in0 = this->mesh.g2l_nodes[this->mesh.faceSub[i].inod_glob[0]];
			in1 = this->mesh.g2l_nodes[this->mesh.faceSub[i].inod_glob[1]];
			in2 = this->mesh.g2l_nodes[this->mesh.faceSub[i].inod_glob[2]];
			AreaOfEl = fabs(this->mesh.coordinateSub[in1].x -
								 this->mesh.coordinateSub[in0].x) *
								 fabs(this->mesh.coordinateSub[in2].y -
								 this->mesh.coordinateSub[in1].y);
			sArea +=    AreaOfEl;
			for (int j = 0;j<4;j++){
				iDof = 3*(this->mesh.g2l_nodes[this->mesh.faceSub[i].inod_glob[j]])+2;
				f[iDof] += 0.25*fz_total*AreaOfEl;
			}
		}
	}
}
//
void CFem::exportRHS(double *f) {
//
//#ifdef flag_Fortran
//  int	startInd_C_or_F = 1;
//#else 
//  int startInd_C_or_F = 0;
//#endif
//	if (this->domainG->i_load_step == 0 && this->domainG->i_sub_step == 0) {
//		this->bound_cond.neuBC->ind = new int[this->domainG->neqAll];
//		this->bound_cond.neuBC->val = new double[this->domainG->neqAll];
//		this->bound_cond.neuBC->n = this->domainG->neqAll;
//	}
//
//	for (int i = 0; i < this->domainG->neqAll; i++) {
//		this->bound_cond.neuBC->ind[i] = i + startInd_C_or_F;
//		this->bound_cond.neuBC->val[i] = f[i];
//	}
}

double CFem::relativeForceEquilib(double *Ku, double *f, int indBC[], int nDir, int ndof) {

	int k = 0;
	double sum_f = 0., sum_Ku_m_f = 0, c;

	for (int i = 0; i < ndof; i++) {
		if (k < nDir && i == indBC[k]) {
			k++;
		} else {
			c = (Ku[i] - f[i]);
			sum_Ku_m_f += c * c;
			sum_f += f[i] * f[i];
		}
	}

	if (sum_f == 0) {
		sum_f = 1.;
	}

	return sqrt(sum_Ku_m_f / sum_f);

}

void CFem::statisticPrintToDisp(CKSparse *KSparse, CDomain *domainG, double *u,
																double * f, double * Ku){
	//int	globoalIndSub = this->i_domGlobal;
	int	globoalIndSub = this->mesh.i_domGlobalMeshGen;
  //MPI_Comm_rank(comm, &MPIrank);


	int &n = domainG->neqAll;
	double * tmp_abs_u = new double[n];
	//
	//
	for (int i = 0;i<n;i++){
		tmp_abs_u[i] = fabs(u[i]);
	}

	qsort(tmp_abs_u, n, sizeof (double), CLinearAlgebra::compareDouble);
	//
	KSparse->multAx(Ku, u, domainG->neqSub[this->i_domOnClust],0);
	//
//	double norm_f_m_Ku = CFem::relativeForceEquilib(Ku,f,
//		this->bound_cond.dirBC->ind,this->bound_cond.dirBC->n,this->domainG->neqAll);
  double norm_f_m_Ku = 1.0e5;
	//
	if (!globoalIndSub) {
		printf("||Ku-f||/||f||=%3.3e\n", norm_f_m_Ku);
		printf("+----%d load step----+\n+----%d  sub step----+\n",
					 domainG->i_load_step, domainG->i_sub_step);
		printf("+              max(|u|)=%3.3e\n",tmp_abs_u[n-1]);
	};
	//
	delete [] tmp_abs_u;
	//
}

void CFem::get_DP_DOFs(CDomain *domainG){
  int tmp_i;
  int I,J,K;
  int nx    = domainG->nx; 
  int ny    = domainG->ny; 
  int nz    = domainG->nz;
  int nxSub = domainG->nxSub; 
  int nySub = domainG->nySub; 
  int nzSub = domainG->nzSub;
  int Nx    = domainG->Nx; 
  int Ny    = domainG->Ny; 
  int Nz    = domainG->Nz;
  int NxClst= domainG->NxClst; 
  int NyClst= domainG->NyClst; 
  int NzClst= domainG->NzClst;
  int remI,remJ,remK;

  bool flag_DP_real_corners=true;

  this->mesh.DP_DOFsCorners.reserve(24);
//
  int cnt=0;
  for (int i=0;i<8;i++){
    tmp_i = this->mesh.l2g_nodes[this->mesh.fixingNodes[i]];
//    if (tmp_i%domainG)
    K = ceil(double (tmp_i+1)/((nx+1)*(ny+1)));
    int inXYplane = (tmp_i+1)-(K-1)*(nx+1)*(ny+1);
    J = ceil(double( inXYplane)/(nx+1));
    I = inXYplane - (nx+1)*(J-1);
    I--; J--; K--;
    remI=I%(nxSub*NxClst);
    remJ=J%(nySub*NyClst);
    remK=K%(nzSub*NzClst);
    if (flag_DP_real_corners && (remI!=0 || remJ!=0 || remK!=0)){
//      printf("isub = %d, tmp = %d,  I = %d, J = %d, K = %d\n",
//             this->mesh.i_domOnClust,tmp_i,I,J,K);  
      this->mesh.DP_DOFsCorners.push_back(3*tmp_i+0);
      this->mesh.DP_DOFsCorners.push_back(3*tmp_i+1);
      this->mesh.DP_DOFsCorners.push_back(3*tmp_i+2);
      cnt+=3;
    }
  }

  int n_DP_DOFsAll=this->mesh.DP_DOFsCorners.size()+3*this->mesh.DP_NodesAll.size();
  this->mesh.DP_DOFsAll.reserve(n_DP_DOFsAll);
  this->mesh.DP_DOFsAll.insert( this->mesh.DP_DOFsAll.begin(),
                                this->mesh.DP_DOFsCorners.begin(),
                                this->mesh.DP_DOFsCorners.end());
  cnt=this->mesh.DP_DOFsCorners.size();
  std::vector<longInt>::iterator it;
  for ( it = this->mesh.DP_NodesAll.begin()  ;
        it != this->mesh.DP_NodesAll.end();it++){
    this->mesh.DP_DOFsAll.push_back(*it*3+0);
    this->mesh.DP_DOFsAll.push_back(*it*3+1);
    this->mesh.DP_DOFsAll.push_back(*it*3+2);
    cnt+=3;
  }

  if (this->mesh.DP_DOFsAll.size() > 0) {
	  qsort(&(this->mesh.DP_DOFsAll[0]),this->mesh.DP_DOFsAll.size(), 
						  sizeof (longInt), CLinearAlgebra::compareLongInt); // POZOR - longInt
	  it = unique (this->mesh.DP_DOFsAll.begin(), this->mesh.DP_DOFsAll.end());
	  this->mesh.DP_DOFsAll.resize( distance(this->mesh.DP_DOFsAll.begin(),it));
  }
//  for ( vector<int>::iterator it = this->mesh.DP_DOFsAll.begin()  ;
//        it != this->mesh.DP_DOFsAll.end();it++)
//        printf("DP_DOFsAll = %d\n",*it);

}
