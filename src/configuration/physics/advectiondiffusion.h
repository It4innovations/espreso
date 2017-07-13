
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_

#include "../material/coordinatesystem.h"
#include "../material/holder.h"
#include "solver.h"
#include "../solver.h"
#include "nonlinearsolver.h"
#include "advectiondiffusionconvection.h"
#include "advectiondiffusionradiation.h"

namespace espreso {

struct AdvectionDiffusionNonLinearConvergence: public NonLinearConvergence {

	virtual bool checkSolution() const { return temperature; }
	virtual bool checkResidual() const { return heat; }

	virtual double requestedSolution() const { return temperature_residual; }
	virtual double requestedResidual() const { return heat_residual; }

	PARAMETER(bool, temperature, "Turn on/off temperature residual check.", true);
	PARAMETER(bool, heat       , "Turn on/off heat residual check."       , false);

	PARAMETER(double, temperature_residual, "Requested temperature residual", 1e-3);
	PARAMETER(double, heat_residual       , "Requested heat residual"       , 1e-3);
};


struct AdvectionDiffusionConfiguration: public Configuration {

	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	OPTION(STABILIZATION, stabilization, "The type of the stabilization.", STABILIZATION::SUPG, OPTIONS({
		{ "SUPG", STABILIZATION::SUPG, "SUPG stabilization." },
		{ "CAU" , STABILIZATION::CAU , "CAU stabilization." },
	}));
	PARAMETER(double, sigma, "Inconsistent stabilization parameters.", 0);
	PARAMETER(bool, tangent_matrix_correction, "Add derivation matrix to stiffness matrix.", 0);

	SUBCONFIG(PhysicsSolver<AdvectionDiffusionNonLinearConvergence>, physics_solver, "Settings of physics solver.");

	SUBMAPTOMAP(size_t, std::string, std::string, heat_flux, "Heat flux", "TIME_STEP", "Heat flux settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, heat_flow, "Heat flow", "TIME_STEP", "Heat flow settings for the load step", "REGION", "EXPRESSION");

	SUBMAPTOMAPTOCONFIG(size_t, std::string, AdvectionDiffusionConvection, convection, "Region with convective heat flux",
			"REGION", "Convection for a given region", "1", "Settings for the load step");
	SUBMAPTOMAPTOCONFIG(size_t, std::string, AdvectionDiffusionRadiation, diffuse_radiation, "Region with diffuse radiation",
			"REGION", "Diffuse radiation for a given region", "1", "Settings for the load step");

	SUBMAP(std::string, std::string, initial_temperature , "Regions initial temperature", "REGION", "EXPRESSION");

	SUBMAPTOMAP(size_t, std::string, std::string, temperature        , "Temperature"       , "TIME_STEP", "Temperature settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, heat_source        , "Heat source"       , "TIME_STEP", "Heat source settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, translation_motions, "Translation motion", "TIME_STEP", "Translation motion settings for the load step", "REGION", "EXPRESSION");

	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "REGION", "MATERIAL");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);

	PARAMETER(bool, bem4i, "Assemble matrices using BEM4I library.", false);
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSION_H_ */
