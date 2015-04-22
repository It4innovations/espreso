#ifndef MESH_H_1
#define MESH_H_1

#include "utility.h"
#include "Solid45.h"
#include "Coordinate.h"
#include "RectangularFace.h"
#include "Domain_g.h"
//#include "Solid45NodNumbLoc.h"

class CMesh {
public:
	  CSolid45 *element; /*G: all elements */
	  CCoordinate *coordinateSub; /*G: coordinates */
	  CRectangularFace *faceSub ; /*G: */
		std::map<longInt,int> g2l; // mapping from global to local numbering
		std::vector<longInt>  l2g; // mapping from local to global numbering
		std::map<longInt,int> g2l_nodes; // mapping from global to local numbering
		std::vector<longInt>  l2g_nodes; // mapping from local to global numbering
		std::vector<longInt>  externalDOFs; // mapping from local to global numbering
		std::vector<longInt>  DP_DOFsCorners; // mapping from local to global numbering
		std::vector<longInt>  DP_NodesAll; // mapping from local to global numbering
		std::vector<longInt>  DP_DOFsAll; // mapping from local to global numbering
    int * nodOnEdg;
    int * neighbSub;
    int * fixingDOFs;
    int * fixingNodes;
    longInt * DP_Nodes;
    int i_domOnClust;
    int i_domGlobal;
    int i_domGlobalMeshGen;
    CDomain *domainG;

public:
    CMesh();
//	  CMesh(CDomain *domainG,int i_domOnClust);
    CMesh(const CMesh& mesh):
        g2l(mesh.g2l), 
        l2g(mesh.l2g), 
        DP_DOFsCorners(mesh.DP_DOFsCorners),
        DP_NodesAll(mesh.DP_NodesAll),
        DP_DOFsAll(mesh.DP_DOFsAll),
        g2l_nodes(mesh.g2l_nodes),
        l2g_nodes(mesh.l2g_nodes),
        externalDOFs(mesh.externalDOFs)
    {
      printf(" ---------------- ---------------- COPY -----------\n");
      if (mesh.domainG!=NULL) {
        domainG = mesh.domainG;
        i_domOnClust = mesh.i_domOnClust;
        i_domGlobal= mesh.i_domGlobal;
        i_domGlobalMeshGen= mesh.i_domGlobalMeshGen;
        int n_nodOnEdg =	
            (domainG->nxSub+1)* (domainG->nySub+1)* (domainG->nzSub+1) - 
            (domainG->nxSub-1)* (domainG->nySub-1)* (domainG->nzSub-1);
        nodOnEdg = new int[n_nodOnEdg];
        memcpy(nodOnEdg,mesh.nodOnEdg,n_nodOnEdg*sizeof(int));
        neighbSub = new int[26];
        memcpy(neighbSub,mesh.neighbSub,26*sizeof(int));
        fixingDOFs = new int[24];
        memcpy(fixingDOFs,mesh.fixingDOFs,24*sizeof(int));
        fixingNodes = new int[8];
        memcpy(fixingNodes,mesh.fixingNodes,8*sizeof(int));
        int tmpInt=domainG->n_nodsSub[i_domOnClust];
        coordinateSub = new CCoordinate[tmpInt];
        memcpy(coordinateSub,mesh.coordinateSub,tmpInt*sizeof(CCoordinate));
        tmpInt=domainG->n_facesSub[i_domOnClust];
        faceSub = new CRectangularFace[tmpInt];
        memcpy(faceSub,mesh.faceSub,tmpInt*sizeof(CRectangularFace));
        tmpInt=domainG->n_elementsSub[i_domOnClust];
        element = new CSolid45[tmpInt];
        memcpy(element,mesh.element,tmpInt*sizeof(CSolid45));
      }
      else{
        domainG=NULL;
        i_domOnClust=0;
        i_domGlobal=0;
        i_domGlobalMeshGen=0;
        element=NULL;
        coordinateSub=NULL;
        faceSub=NULL;
        nodOnEdg=NULL;
        neighbSub=NULL;
        fixingDOFs=NULL;
        fixingNodes=NULL;

      }
    }

    CMesh& operator=( const CMesh& mesh) {
// ----
 
    {
      if (this==&mesh) {return *this;}
      g2l=mesh.g2l; 
      l2g=mesh.l2g; 
      DP_DOFsCorners=mesh.DP_DOFsCorners;
      DP_NodesAll=mesh.DP_NodesAll;
      DP_DOFsAll=mesh.DP_DOFsAll;
      g2l_nodes=mesh.g2l_nodes;
      l2g_nodes=mesh.l2g_nodes;
      externalDOFs=mesh.externalDOFs;
      printf(" ---------------- ---------------- ASSIG-----------\n");
// ---------------------------------------
      if (mesh.domainG!=NULL) {
        domainG = mesh.domainG;
        i_domOnClust = mesh.i_domOnClust;
        i_domGlobal= mesh.i_domGlobal;
        i_domGlobalMeshGen= mesh.i_domGlobalMeshGen;
        int n_nodOnEdg =	
            (domainG->nxSub+1)* (domainG->nySub+1)* (domainG->nzSub+1) - 
            (domainG->nxSub-1)* (domainG->nySub-1)* (domainG->nzSub-1);
        nodOnEdg = new int[n_nodOnEdg];
        memcpy(nodOnEdg,mesh.nodOnEdg,n_nodOnEdg*sizeof(int));
        neighbSub = new int[26];
        memcpy(neighbSub,mesh.neighbSub,26*sizeof(int));
        fixingDOFs = new int[24];
        memcpy(fixingDOFs,mesh.fixingDOFs,24*sizeof(int));
        fixingNodes = new int[8];
        memcpy(fixingNodes,mesh.fixingNodes,8*sizeof(int));
        int tmpInt=domainG->n_nodsSub[i_domOnClust];
        coordinateSub = new CCoordinate[tmpInt];
        memcpy(coordinateSub,mesh.coordinateSub,tmpInt*sizeof(CCoordinate));
        tmpInt=domainG->n_facesSub[i_domOnClust];
        faceSub = new CRectangularFace[tmpInt];
        memcpy(faceSub,mesh.faceSub,tmpInt*sizeof(CRectangularFace));
        tmpInt=domainG->n_elementsSub[i_domOnClust];
        element = new CSolid45[tmpInt];
        memcpy(element,mesh.element,tmpInt*sizeof(CSolid45));
      }
      else{
        domainG=NULL;
        i_domOnClust=0;
        i_domGlobal=0;
        i_domGlobalMeshGen=0;
        element=NULL;
        coordinateSub=NULL;
        faceSub=NULL;
        nodOnEdg=NULL;
        neighbSub=NULL;
        fixingDOFs=NULL;
        fixingNodes=NULL;

      }
// ---------------------------------------
    }
    }
	virtual ~CMesh();

	void initialize(CDomain *domainG, int i_domOnClust, int i_domGlobal, int i_domGlobalMeshGen);
};

#endif /* MESH_H_1 */
