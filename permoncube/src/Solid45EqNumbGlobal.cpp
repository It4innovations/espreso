#include "Solid45EqNumbGlobal.h"

CSolid45EqNumbGlobal::CSolid45EqNumbGlobal()
{
  memset(ieq,0, 24 * sizeof(longInt));
  stif_loc = new CStiffnessLocal();
}

CSolid45EqNumbGlobal::~CSolid45EqNumbGlobal() {
  if (stif_loc!=NULL){ 
    delete stif_loc;
    stif_loc=NULL;
  }
}


