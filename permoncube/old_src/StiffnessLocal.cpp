
#include "StiffnessLocal.h"

CStiffnessLocal::CStiffnessLocal()
{
  //double value_K[576]; /*{x1,x2,x3, ...x8,y1,y2,y3,.....z6,z7,z8}*/
  //	
  value_K = new double[576];
  memset(value_K,0, 576 * sizeof(double));
  value_f = new double[24];
  memset(value_f,0, 24 * sizeof(double));
  value_f_internal = new double[24];
  memset(value_f_internal,0, 24 * sizeof(double));
  value_u = new double[24];
  memset(value_u,0, 24 * sizeof(double));
  value_du = new double[24];
  memset(value_du,0, 24 * sizeof(double));
  // *** plasticity part ***
  if (!gv_flag_linear_system) {
    currentKinematicHardening = new double[6*8];
    memset(currentKinematicHardening,0, 6*8 * sizeof(double));

    deltaCurrentIsotropicHardening = new double[8];
    memset(deltaCurrentIsotropicHardening, 0,8 * sizeof(double));

    deltaCurrentKinematicHardening = new double[6*8];
    memset(deltaCurrentKinematicHardening,0, 6*8 * sizeof(double));

    sigepl = new double[8];
    memset(sigepl, 0,8 * sizeof(double));

    sigrat = new double[8];
    memset(sigrat, 0,8 * sizeof(double));

    depeq = new double[8];
    memset(depeq, 0,8 * sizeof(double));

    stress = new double[6*8];
    memset(stress, 0,6*8 * sizeof(double));

    tmp_stress = new double[6*8];
    memset(tmp_stress, 0,6*8 * sizeof(double));

    epel = new double[6*8];
    memset(epel, 0,6*8 * sizeof(double));   

    tmp_epel = new double[6*8];
    memset(tmp_epel, 0,6*8 * sizeof(double));   

    eppl = new double[6*8];
    memset(eppl, 0,6*8 * sizeof(double));

    tmp_eppl = new double[6*8];
    memset(tmp_eppl, 0,6*8 * sizeof(double));

    statev = new double[6*8*8];
    memset(statev, 0,6*6*8 * sizeof(double));

    tmp_statev = new double[6*8*8];
    memset(tmp_statev, 0,6*6*8 * sizeof(double));

    epeq = new double[8];
    memset(epeq, 0,8 * sizeof(double));

    tmp_epeq = new double[8];
    memset(tmp_epeq, 0,8 * sizeof(double));

    plwork = new double[8];
    memset(plwork, 0,8 * sizeof(double));

    tmp_plwork = new double[8];
    memset(tmp_plwork, 0,8 * sizeof(double));

    currentIsotropicHardening = new double[8];
    tmp_currentIsotropicHardening = new double[8];
    for (int i = 0;i<8;i++){
      currentIsotropicHardening[i] = 450.;
      tmp_currentIsotropicHardening[i] = 450.;

    }
  }
  else {
    currentKinematicHardening = NULL;
    deltaCurrentIsotropicHardening = NULL;
    deltaCurrentKinematicHardening = NULL;
    sigepl = NULL;
    sigrat = NULL;
    depeq = NULL;
    stress = NULL;
    tmp_stress = NULL;
    currentIsotropicHardening = NULL;
    tmp_currentIsotropicHardening = NULL;
    epel = NULL;
    tmp_epel  = NULL;
    eppl = NULL;
    tmp_eppl = NULL;
    statev = NULL;
    tmp_statev = NULL;
    epeq = NULL;
    tmp_epeq = NULL;
    plwork = NULL;
    tmp_plwork = NULL;
  }
}

CStiffnessLocal::~CStiffnessLocal() {
  if (value_K) 											  {	delete [] value_K;                          value_K=NULL;} 											  
  if (value_f)  											{	delete [] value_f;                          value_f=NULL;}  											
  if (value_f_internal)						  	{	delete [] value_f_internal;                 value_f_internal=NULL;}						  	
  if (value_u)  											{	delete [] value_u;                          value_u=NULL;}  											
  if (value_du) 											{	delete [] value_du;                         value_du=NULL ;}											
  //                                                                                                                            
  if (currentKinematicHardening) 	  	{	delete [] currentKinematicHardening;        currentKinematicHardening=NULL;} 	  	
  if (deltaCurrentIsotropicHardening) {	delete [] deltaCurrentIsotropicHardening;   deltaCurrentIsotropicHardening=NULL;} 
  if (deltaCurrentKinematicHardening) { delete [] deltaCurrentKinematicHardening;   deltaCurrentKinematicHardening=NULL ;}
  if (sigepl) 												{	delete [] sigepl;                           sigepl=NULL;} 												
  if (sigrat) 												{	delete [] sigrat;                           sigrat=NULL;} 												
  if (depeq) 												  {	delete [] depeq;                            depeq=NULL;} 												  
  if (stress) 												{	delete [] stress;                           stress=NULL ;}												
  if (tmp_stress) 										{	delete [] tmp_stress;                       tmp_stress=NULL;} 										
  if (currentIsotropicHardening) 		  {	delete [] currentIsotropicHardening;        currentIsotropicHardening=NULL;} 		  
  if (tmp_currentIsotropicHardening)  {	delete [] tmp_currentIsotropicHardening;    tmp_currentIsotropicHardening=NULL;}  
  if (epel) 													{	delete [] epel;                             epel=NULL;} 													
  if (tmp_epel) 											{	delete [] tmp_epel;                         tmp_epel=NULL;} 											
  if (eppl) 													{	delete [] eppl;                             eppl=NULL;} 													
  if (tmp_eppl) 											{	delete [] tmp_eppl;                         tmp_eppl=NULL;} 											
  if (statev) 												{	delete [] statev;                           statev=NULL;} 												
  if (tmp_statev) 										{	delete [] tmp_statev;                       tmp_statev=NULL;} 										
  if (epeq) 													{	delete [] epeq;                             epeq=NULL;} 													
  if (tmp_epeq) 											{	delete [] tmp_epeq;                         tmp_epeq=NULL ;}											
  if (plwork) 												{	delete [] plwork;                           plwork=NULL;} 												
  if (tmp_plwork) 										{	delete [] tmp_plwork;                       tmp_plwork=NULL ;}										
}

