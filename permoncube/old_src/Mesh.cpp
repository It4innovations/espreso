#include "Mesh.h"
//
//CMesh::CMesh(CDomain *_domainG,int _i_domOnClust) :
CMesh::CMesh() :
  domainG(NULL),
  i_domOnClust(NULL),
  i_domGlobal(NULL),
  i_domGlobalMeshGen(NULL),
  element(NULL),
  coordinateSub(NULL),
  faceSub(NULL), 
  nodOnEdg(NULL),
  neighbSub(NULL),
  fixingDOFs(NULL)
{}
//
CMesh::~CMesh() {
  if (element)         { delete [] element;        element=NULL; }
  if (coordinateSub)   { delete [] coordinateSub;  coordinateSub=NULL;}
  if (faceSub)         { delete [] faceSub;        faceSub=NULL; }
  if (nodOnEdg)        { delete [] nodOnEdg;       nodOnEdg=NULL;}
  if (neighbSub)       { delete [] neighbSub;      neighbSub=NULL;}
  if (fixingDOFs)      { delete [] fixingDOFs;     fixingDOFs=NULL;}
}
//
void CMesh::initialize(CDomain *domainG,int i_domOnClust,int i_domGlobal,int i_domGlobalMeshGen)
{          
  this->domainG = domainG;
  this->i_domOnClust=i_domOnClust;
  this->i_domGlobal=i_domGlobal;
  this->i_domGlobalMeshGen=i_domGlobalMeshGen;
  /*on all nodes*/
  element 			= new CSolid45[domainG->n_elementsSub[i_domOnClust]];
  coordinateSub = new CCoordinate[domainG->n_nodsSub[i_domOnClust]];     
  faceSub 			= new CRectangularFace[domainG->n_facesSub[i_domOnClust]];
  int n_nodOnEdg =	
    (domainG->nxSub+1)* (domainG->nySub+1)* (domainG->nzSub+1) - 
    (domainG->nxSub-1)* (domainG->nySub-1)* (domainG->nzSub-1);
  domainG->n_nodOnEdg[i_domOnClust] = n_nodOnEdg;
  domainG->n_exterDOFs[i_domOnClust] = n_nodOnEdg*3;
  nodOnEdg = new int[n_nodOnEdg];
}
//
