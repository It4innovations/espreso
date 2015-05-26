
#ifndef STIFFNESSLOCAL_H_
#define STIFFNESSLOCAL_H_

#include "utility.h"
extern bool gv_flag_linear_system;

class CStiffnessLocal {
public:
//  {x1,x2,x3, ...x8,y1,y2,y3,.....z6,z7,z8}
  double *value_K; 
  double *value_f; 
  double *value_f_internal; 
  double *value_u; 
  double *value_du;
// *** plasticity part ***
  double *currentKinematicHardening;
  double *deltaCurrentIsotropicHardening; 
  double *deltaCurrentKinematicHardening;
  double *sigepl;
  double *sigrat;
  double *depeq;
  double *stress; 
  double *tmp_stress; 
  double *currentIsotropicHardening; 
  double *tmp_currentIsotropicHardening; 
  double *epel; 
  double *tmp_epel; 
  double *eppl; 
  double *tmp_eppl; 
  double *statev; 
  double *tmp_statev; 
  double *epeq; 
  double *tmp_epeq; 
  double *plwork; 
  double *tmp_plwork; 
    
public:
	  CStiffnessLocal();
	virtual ~CStiffnessLocal();
};

#endif /* STIFFNESSLOCAL_H_ */
