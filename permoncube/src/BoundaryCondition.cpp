
#include "BoundaryCondition.h"

CBoundaryCondition::CBoundaryCondition()
{
	dirBCSub	    = new CDirichletBC();
	dirBC_global  = new CDirichletBC();
	neuBC			    = new CNeumannBC();
	conBCSub	    = new CContactBC();
}

CBoundaryCondition::CBoundaryCondition (const CBoundaryCondition & _BoundaryCondition)
{
  dirBCSub = new CDirichletBC(*_BoundaryCondition.dirBCSub);
  // dopsat kopirovaci konstruktor dirBCSub

}

CBoundaryCondition::~CBoundaryCondition()
{
  //printf(" ++++ ++++ +++ ++++ ++++ in destructor of BoundaryCondition\n");
	delete dirBCSub;
  if (dirBC_global) { delete dirBC_global; dirBC_global=NULL; }
	if (neuBC)        { delete neuBC;        neuBC=NULL;}
	if (conBCSub)     { delete conBCSub;     conBCSub=NULL;}
}

