
#include "ElementLibrary.h"

CElementLibrary::CElementLibrary() {
}

CElementLibrary::~CElementLibrary() {
}

void CElementLibrary::modulus(int i_sub_step,
        double *c,double E, double mu, 
				double *epel, 
				double *eppl,
        double *currentIsotropicHardening,
        double *deltaCurrentIsotropicHardening,
        double *currentKinematicHardening,
        double *deltaCurrentKinematicHardening,
        double *Gt, 
        double *u,double *du, double *stress, 
        double *ct)
{
// TODO material parameters should be defined by user-------------------
	double Hm = 1.e5;
//double sY = 450.;
// ---------------------------------------------------------------------
	double g   = 0.5*E/(1.+mu); // Lamme constant (g = mu)
  double con2;
  //TODO depel is not obligate entity
  // better to pass epel = epel_prew + depel
  double depel[6];
  memset(depel,0, 6 * sizeof(double));


  for (int i = 0;i<6;i++) {
    for (int j = 0;j<24;j++) {
      depel[i] += Gt[i*24+j]*du[j];
    } 
    epel[i]+=depel[i];
  } 

  double ep_plus_dep[6];
  memset(ep_plus_dep,0, 6 * sizeof(double));
  for (int i = 0;i<6;i++) {
    for (int j = 0;j<24;j++) {
      ep_plus_dep[i] += Gt[i*24+j]*(u[j] + du[j]);
    } 
  } 
  
  double sigma_c[6];  
  memset(sigma_c,0, 6 * sizeof(double));
  for (int i = 0;i<6;i++) {
    for (int j = 0;j<6;j++) {
      sigma_c[i] += c[i*6+j]*ep_plus_dep[j];
    } 
  } 
  
  double sigma_p[6];  
  memset(sigma_p,0, 6 * sizeof(double));
  for (int i = 0;i<6;i++) {
    for (int j = 0;j<6;j++) {
      sigma_p[i] += c[i*6+j]*eppl[j];
    } 
  }

  double sigma_t[6];  
  for (int i = 0;i<6;i++) {
    sigma_t[i] = sigma_c[i] - sigma_p[i];
  } 


  double dev_sigma_t[6];
  dev_sigma_t[0] = (sigma_t[0]*2-sigma_t[1]-sigma_t[2])/3.;
  dev_sigma_t[1] = (-sigma_t[0]+sigma_t[1]*2-sigma_t[2])/3.;
  dev_sigma_t[2] = (-sigma_t[0]-sigma_t[1]+sigma_t[2]*2)/3.;
  dev_sigma_t[3] = sigma_t[3];
  dev_sigma_t[4] = sigma_t[4];
  dev_sigma_t[5] = sigma_t[5];


  double norm_dev_sigma_t = sqrt(
					  	pow(dev_sigma_t[0],2)+
					  	pow(dev_sigma_t[1],2)+
				  		pow(dev_sigma_t[2],2)+ 
						2*pow(dev_sigma_t[3],2)+ 
						2*pow(dev_sigma_t[4],2)+ 
						2*pow(dev_sigma_t[5],2));

  double hardening_local;
  hardening_local = currentIsotropicHardening[0];
  double fun_plas = sqrt(3./2.)*norm_dev_sigma_t-hardening_local;
  double eps_tau = 1.e-4;

	if (fun_plas <= -eps_tau){
    memcpy(stress,sigma_t,6*sizeof(double));
    memcpy(ct,c,36*sizeof(double));
	}
	else {
    double conIsotrop1 = 3.*g/(3.*g+Hm);
    double n_dev[6], alpha1[6],deppl[6]; 
    double con3 = sqrt(2./3.);
    for (int i = 0;i<6;i++) {
      n_dev[i] = dev_sigma_t[i]/norm_dev_sigma_t;
      alpha1[i] = con3*fun_plas*n_dev[i];
      stress[i] = sigma_c[i] - sigma_p[i] - conIsotrop1*alpha1[i];
      deppl[i] = (1./(2.*g))*conIsotrop1*alpha1[i];
      epel[i] -= deppl[i];
      eppl[i] += deppl[i];
    }

		double kappaLocal = fun_plas/(3.*g+Hm);
		currentIsotropicHardening[0] += kappaLocal*Hm;
     
    double N[36];
    for (int i = 0;i<6;i++) {
		  for (int j = 0;j<6;j++) {
        N[i*6+j] = n_dev[i]*n_dev[j];  
      }
    }
    con2 = sqrt(2./3.)*(hardening_local/norm_dev_sigma_t);
    
		double cd1 = 2./3.,cd2=-1./3.,cd3=0.5;
    double DevMatrixDer[36] =  { cd1,cd2,cd2,0.,0.,0.,
															   cd2,cd1,cd2,0.,0.,0.,
															   cd2,cd2,cd1,0.,0.,0.,
															   0., 0.,0.,cd3 ,0.,0.,
															   0., 0.,0.,0.,cd3 ,0.,
															   0., 0.,0.,0.,0.,cd3};
//
    for (int i = 0;i<6;i++) {
		  for (int j = 0;j<6;j++) {
        ct[i*6+j] = c[i*6+j]
              - 2.*g*conIsotrop1*(DevMatrixDer[i*6+j]
              - con2*(DevMatrixDer[i*6+j]-N[i*6+j]));
      }
    }
  }
}



void CElementLibrary::stima_solid45(
    int i_sub_step,
    double e, double mu, 
		CStiffnessLocal *stif_loc,
		CCoordinate *coordinate, const double *R, double *C,
		double *acceleration, double density)
{
	/* hexahedron
	 * local element matrix
	 * hexahedron (in ANSYS - solid 45)
	 * 8 nodes,3 degrees of freedom per node
	 * 24x24, 576 entries
	 */
	//double density = 2.;
	double r, s, t;
	double Gt[6 * 24];
	double N[8], dNr[8], dNs[8], dNt[8]; // dNx[8], dNy[8], dNz[8]
	double J[9], iJ[9];
	double CG_tmp[6 * 24];
	double dNx_j, dNy_j, dNz_j;
  double ct[36];

	memset(stif_loc->value_f,0,24*sizeof(double));
	memset(stif_loc->value_f_internal,0,24*sizeof(double));
	memset(stif_loc->value_K,0,576*sizeof(double));
	memset(Gt,0,144*sizeof(double));
	memset(N,0,8*sizeof(double));


  // copy of variables over element +++++++++++++++++++++++++++++++++++++++++++
  if (!gv_flag_linear_system){
    memcpy(&stif_loc->tmp_epel[0],&stif_loc->epel[0],48*sizeof(double));
    memcpy(&stif_loc->tmp_eppl[0],&stif_loc->eppl[0],48*sizeof(double));
    memcpy(&stif_loc->tmp_stress[0],&stif_loc->stress[0],48*sizeof(double));
    memcpy(&stif_loc->tmp_currentIsotropicHardening[0],
           &stif_loc->currentIsotropicHardening[0],8*sizeof(double));
    memcpy(&stif_loc->tmp_statev[0],&stif_loc->statev[0],8*6*6*sizeof(double));
    // new ones
    memcpy(&stif_loc->tmp_epeq[0],&stif_loc->epeq[0],8*sizeof(double));
    memcpy(&stif_loc->tmp_plwork[0],&stif_loc->plwork[0],8*sizeof(double));
  }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	for (int i = 0; i < 8; i++) {
		//
		r = R[i];
		s = R[i + 8];
		t = R[i + 16];
		//
		dNr[0] = 0.125 * -(1. - s) * (1. - t);
		dNr[1] = 0.125 * (1. - s) * (1. - t);
		dNr[2] = 0.125 * (1. + s) * (1. - t);
		dNr[3] = 0.125 * -(1. + s) * (1. - t);
		dNr[4] = 0.125 * -(1. - s) * (1. + t);
		dNr[5] = 0.125 * (1. - s) * (1. + t);
		dNr[6] = 0.125 * (1. + s) * (1. + t);
		dNr[7] = 0.125 * -(1. + s) * (1. + t);
		//
		dNs[0] = 0.125 * -(1. - r) * (1. - t);
		dNs[1] = 0.125 * -(1. + r) * (1. - t);
		dNs[2] = 0.125 * (1. + r) * (1. - t);
		dNs[3] = 0.125 * (1. - r) * (1. - t);
		dNs[4] = 0.125 * -(1. - r) * (1. + t);
		dNs[5] = 0.125 * -(1. + r) * (1. + t);
		dNs[6] = 0.125 * (1. + r) * (1. + t);
		dNs[7] = 0.125 * (1. - r) * (1. + t);
		//
		dNt[0] = 0.125 * -(1. - r) * (1. - s);
		dNt[1] = 0.125 * -(1. + r) * (1. - s);
		dNt[2] = 0.125 * -(1. + r) * (1. + s);
		dNt[3] = 0.125 * -(1. - r) * (1. + s);
		dNt[4] = 0.125 * (1. - r) * (1. - s);
		dNt[5] = 0.125 * (1. + r) * (1. - s);
		dNt[6] = 0.125 * (1. + r) * (1. + s);
		dNt[7] = 0.125 * (1. - r) * (1. + s);
		//
		N[0] = 0.125 * (1. - r) * (1. - s) * (1. - t);
		N[1] = 0.125 * (r + 1.) * (1. - s) * (1. - t);
		N[2] = 0.125 * (r + 1.) * (s + 1.) * (1. - t);
		N[3] = 0.125 * (1. - r) * (s + 1.) * (1. - t);
		N[4] = 0.125 * (1. - r) * (1. - s) * (t + 1.);
		N[5] = 0.125 * (r + 1.) * (1. - s) * (t + 1.);
		N[6] = 0.125 * (r + 1.) * (s + 1.) * (t + 1.);
		N[7] = 0.125 * (1. - r) * (s + 1.) * (t + 1.);
		//
		for (int j = 0; j < 9; j++) {
			J[j] = 0.;
			iJ[j] = 0.;
		}
		//
		for (int j = 0; j < 8; j++) {
			J[0] += (dNr[j] * coordinate[j].x);
			J[1] += (dNs[j] * coordinate[j].x);
			J[2] += (dNt[j] * coordinate[j].x);
			J[3] += (dNr[j] * coordinate[j].y);
			J[4] += (dNs[j] * coordinate[j].y);
			J[5] += (dNt[j] * coordinate[j].y);
			J[6] += (dNr[j] * coordinate[j].z);
			J[7] += (dNs[j] * coordinate[j].z);
			J[8] += (dNt[j] * coordinate[j].z);
		}
		//
		double detJ = CLinearAlgebra::inverse_matrix_3x3(J, iJ);
		//
		for (int j = 0; j < 8; j++) {
			dNx_j = iJ[0] * dNr[j] + iJ[3] * dNs[j] + iJ[6] * dNt[j];
			dNy_j = iJ[1] * dNr[j] + iJ[4] * dNs[j] + iJ[7] * dNt[j];
			dNz_j = iJ[2] * dNr[j] + iJ[5] * dNs[j] + iJ[8] * dNt[j];
			Gt[j] = dNx_j; // x
			Gt[j + 32] = dNy_j; // y
			Gt[j + 64] = dNz_j; // z
			Gt[j + 72] = dNy_j; // y-x
			Gt[j + 80] = dNx_j; // y-x
			Gt[j + 104] = dNz_j; // y-z
			Gt[j + 112] = dNy_j; // y-z
			Gt[j + 120] = dNz_j; // x-z
			Gt[j + 136] = dNx_j; // x-z
		}
		/*
		for (int j = 0; j < (6 * 24); j++) {
			CG_tmp[j] = 0.;
		}
    */    
	  memset(CG_tmp,0,144*sizeof(double));
    //
    if (gv_flag_linear_system){
      memcpy(ct,C,36*sizeof(double));
    }
    else {
  //#define usermat_fortran
#ifdef usermat_fortran
    int elem=1,mat=1,ncomp=6,kfirst=0,kfsteq=0,ktform=1;
    double *prop, *tmpdouble;
    double timval,timinc,tem,dtem,toffst,flu,dflu;
 
    /* 1 - has to be passed
     * 0 - does not have to be */
//    double ep_plus_dep[6];

    for (int I = 0;I<6;I++) {
    //  stif_loc->tmp_epel[i*6+I]=0.0;
      for (int J = 0;J<24;J++) {
        stif_loc->tmp_epel[i*6+I] +=
            Gt[I*24+J]*(stif_loc->value_u[J] + stif_loc->value_du[J]);
        } 
      } 
    userpl_(elem,                                  // 0  
            i,                                     // 0 
            mat,                                   // 0
            ncomp,                                 // 1
            kfirst,                                // 0
            kfsteq,                                // 0
            e,                                     // 1
            mu,                                    // 1
            density,                               // 0
            prop,                                  // 0
            C,                                     // 1
            ktform,                                //-1 
            timval,                                // 0
            timinc,                                // 0
            tem,                                   // 0
            dtem,                                  // 0
            toffst,                                // 0
            flu,                                   // 0
            dflu,                                  // 0
            &(stif_loc->tmp_epel[i*6]),            // 1 INOUT
            &(stif_loc->tmp_eppl[i*6]),            // 1 INOUT
            &(stif_loc->tmp_statev[i*6]),          // 1 INOUT
            tmpdouble,//&(stif_loc->tmp_usvr[i]),  // 0 INOUT
            (stif_loc->tmp_epeq[i]),               // 1 INOUT
            (stif_loc->tmp_plwork[i]),             // 1 INOUT
            (stif_loc->sigepl[i]),                 // 1 INOUT
            (stif_loc->sigrat[i]),                 // 1 OUT
            (stif_loc->depeq[i]),                  // 1 INOUT
            ct);                                   // 1 INOUT
      
    for (int I=0;I<6;I++){
      stif_loc->tmp_stress[i*6+I] = 0.0;
      for (int J=0;J<6;J++){
        stif_loc->tmp_stress[i*6+I]+=C[I+J*6]*stif_loc->tmp_epel[i*6+J];
      }
      }
#else
      
      CElementLibrary::modulus(i_sub_step, C,e, mu, 
                    &(stif_loc->tmp_epel[i*6]), 
                    &(stif_loc->tmp_eppl[i*6]),		
                    &(stif_loc->tmp_currentIsotropicHardening[i]),
                    &(stif_loc->deltaCurrentIsotropicHardening[i]),
                    &(stif_loc->currentKinematicHardening[i*6]),		
                    &(stif_loc->deltaCurrentKinematicHardening[i*6]),		
                    Gt, 
                    stif_loc->value_u,
                    stif_loc->value_du,
                    &(stif_loc->tmp_stress[i*6]),
                    ct);
#endif
      //
      for (int I=0;I<24;I++){
        for (int J = 0;J<6;J++){
          stif_loc->value_f_internal[I] += 
          Gt[I+J*24]*stif_loc->tmp_stress[i*6+J]* detJ;
        }
      }
    } //if (flag_linear_system){
//
    for (int I = 0; I < 6; I++) {
      for (int J = 0; J < 24; J++) {
        CG_tmp[I * 24 + J] = 0.;
        for (int K = 0; K < 6; K++) {
          CG_tmp[I * 24 + J] += (ct[I + K * 6] * Gt[24 * K + J]);
        }
      }
    }
    for (int I = 0; I < 24; I++) {
      for (int J = 0; J < 24; J++) {
        for (int K = 0; K < 6; K++) {
          stif_loc->value_K[I + J * 24] += (Gt[I + K * 24]
              * CG_tmp[K * 24 + J] * detJ);
        }
      }
    }
    //
    for (int j = 0; j < 8; j++) {
      stif_loc->value_f[j] += acceleration[0] * density * N[j] * detJ;
      stif_loc->value_f[j + 8] += acceleration[1] * density * N[j] * detJ;
      stif_loc->value_f[j + 16] += acceleration[2] * density * N[j]
          * detJ;
    }
	}
}
