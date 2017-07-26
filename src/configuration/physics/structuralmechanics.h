
#ifndef SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS_H_

#include "../material/coordinatesystem.h"
#include "../material/holder.h"
#include "solver.h"
#include "../solver.h"
#include "nonlinearsolver.h"

namespace espreso {

struct StructuralMechanicsNonLinearConvergence: public NonLinearConvergence {

	virtual bool checkSolution() const { return displacement; }
	virtual bool checkResidual() const { return give_me_name; }

	virtual double requestedSolution() const { return displacement_residual; }
	virtual double requestedResidual() const { return give_me_name_residual; }

	PARAMETER(bool, displacement, "Turn on/off displacement residual check.", true);
	PARAMETER(bool, give_me_name, "Turn on/off give_me_name residual check.", false);

	PARAMETER(double, displacement_residual, "Requested displacement residual", 1e-3);
	PARAMETER(double, give_me_name_residual, "Requested give_me_name residual", 1e-3);
};

struct StructuralMechanicsConfiguration: public Configuration {

	SUBCONFIG(PhysicsSolver<StructuralMechanicsNonLinearConvergence>, physics_solver, "Settings of physics solver.");

	SUBMAPTOMAP(size_t, std::string, std::string, displacement       , "Displacement"       , "TIME_STEP", "Displacement settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, normal_presure     , "Normal presure"     , "TIME_STEP", "Normal presure settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, initial_temperature, "Initial temperature", "TIME_STEP", "Initial temperature settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, temperature        , "Temperature"        , "TIME_STEP", "Temperature settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, acceleration       , "Acceleration"       , "TIME_STEP", "Acceleration settings for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, obstacle           , "Obstacle"           , "TIME_STEP", "Obstacles for the load step", "REGION", "EXPRESSION");
	SUBMAPTOMAP(size_t, std::string, std::string, normal_direction   , "Normal Direction"   , "TIME_STEP", "Normal Direction of obstacle settings for the load step", "REGION", "EXPRESSION");

	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "REGION", "MATERIAL");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);

	PARAMETER(bool, bem4i, "Assemble matrices using BEM4I library.", false);
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS_H_ */
