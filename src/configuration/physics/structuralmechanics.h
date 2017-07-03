
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

	SUBMAPTOMAP(size_t, std::string, std::string, displacement       , "Displacement"       , "1", "Displacement settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, normal_presure     , "Normal presure"     , "1", "Normal presure settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, initial_temperature, "Initial temperature", "1", "Initial temperature settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, temperature        , "Temperature"        , "1", "Temperature settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, acceleration       , "Acceleration"       , "1", "Acceleration settings for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, obstacle           , "Obstacle"           , "1", "Obstacles for load step '1'", "<REGION>", "<EXPRESSION>");
	SUBMAPTOMAP(size_t, std::string, std::string, normal_direction   , "Normal Direction"   , "1", "Normal Direction of obstacle settings for load step '1'", "<REGION>", "<EXPRESSION>");

	SUBMAP(std::string, std::string, material_set, "Assign materials to regions", "<REGION>", "<MATERIAL_NAME>");

	PARAMETER(bool, post_process, "Turn on/off results post processing.", true);

	PARAMETER(bool, bem4i, "Assemble matrices using BEM4I library.", false);
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_STRUCTURALMECHANICS_H_ */
