
#ifndef ELEMENTLIBRARY_H_
#define ELEMENTLIBRARY_H_

#include "StiffnessLocal.h"
#include "Coordinate.h"
#include "LinearAlgebra.h"

extern bool gv_flag_linear_system;
class CElementLibrary {
public:
	CElementLibrary();
	virtual ~CElementLibrary();
public:
	static void stima_solid45(int i_sub_step, double E, double mu,
		  CStiffnessLocal *stif_loc,
			CCoordinate *coordinate, const double *R, double *C,
			double *acceleration, double density);
  static void modulus(int i_sub_step, double *C,double E, double mu, 
				double *epel,
				double *eppl,
        double *currentIsotropicHardening,
        double *deltaCurrentIsotropicHardening,
        double *currentKinematicHardening,
        double *deltaCurrentKinematicHardening,
        double *Gt, 
        double *u,double *du, double *stress, 
        double *ct);

};


#endif /* ELEMENTLIBRARY_H_ */
