
#ifndef SRC_CONFIGURATION_GLOBALCONFIGURATION_H_
#define SRC_CONFIGURATION_GLOBALCONFIGURATION_H_

#include "../configuration/decomposer.h"
#include "../configuration/environment.h"
#include "../configuration/output.h"
#include "../configuration/physics.h"
#include "../configuration/results.h"
#include "input/input.h"
#include "input/inputgenerator.h"
#include "material/holder.h"
#include "physics/advectiondiffusion2d.h"
#include "physics/advectiondiffusion3d.h"
#include "physics/structuralmechanics2d.h"
#include "physics/structuralmechanics3d.h"
#include "physics/shallowwater2d.h"
#include "reader/reader.h"

namespace espreso {

enum class INPUT {
	WORKBENCH = 0,
	OPENFOAM = 1,
	ESDATA = 2,
	GENERATOR = 3
};

struct GlobalConfiguration: public Configuration {

	GlobalConfiguration(const std::string &file) { Reader::read(*this, file, this->default_args, this->variables); Reader::set(this->env, this->output); }
	GlobalConfiguration(int *argc, char ***argv) { Reader::read(*this, argc, argv, this->default_args, this->variables); Reader::set(this->env, this->output); }

	void print() { Reader::print(*this); }
	void store() { Reader::store(*this, { ".*" }); }

	SUBMAP(size_t, std::string, default_args, "List of default values for arguments - [ARG*].", "ARG", "EXPRESSION");
	SUBMAP(std::string, std::string, variables, "List of variables usable in *.ecf file.", "VARIABLE", "EXPRESSION");

	SUBCONFIG(Environment        , env         , "Environment dependent variables (set by ./env/threading.* scripts).");
	SUBCONFIG(OutputConfiguration, output      , "Output settings.");

	SUBCONFIG(ESPRESOGenerator   , generator   , "ESPRESO internal mesh generator.");
	SUBCONFIG(ESPRESOInput       , workbench   , "Mesh description in Ansys Workbench format.");
	SUBCONFIG(ESPRESOInput       , openfoam    , "Mesh description in OpenFOAM format.");
	SUBCONFIG(ESPRESOInput       , esdata      , "Mesh description in ESPRESO internal binary format.");

	SUBCONFIG(AdvectionDiffusion2DConfiguration , advection_diffusion_2D , "2D advection diffusiuon solver.");
	SUBCONFIG(AdvectionDiffusion3DConfiguration , advection_diffusion_3D , "3D advection diffusiuon solver.");
	SUBCONFIG(StructuralMechanics2DConfiguration, structural_mechanics_2D, "2D structural mechanics solver.");
	SUBCONFIG(StructuralMechanics3DConfiguration, structural_mechanics_3D, "3D structural mechanics solver.");
	SUBCONFIG(ShallowWater2DConfiguration       , shallow_water_2D       , "2D shallow water solver.");

	SUBCONFIG(Results, results, "Expected output results.");

	SUBCONFIG(Decomposer, decomposer, "./decomposer configuration");

	OPTION(INPUT, input, "test input", INPUT::GENERATOR, OPTIONS({
			{ "WORKBENCH", INPUT::WORKBENCH, { "workbench" }, "Ansys Workbench input file" },
			{ "OPENFOAM", INPUT::OPENFOAM, { "openfoam" }, "OpenFOAM input format" },
			{ "ESDATA", INPUT::ESDATA, { "esdata" }, "ESPRESO binary format" },
			{ "GENERATOR", INPUT::GENERATOR, { "generator" }, "ESPRESO internal generator" }
	}));

	OPTION(PHYSICS, physics, "Used physics", PHYSICS::STRUCTURAL_MECHANICS_3D, OPTIONS({
		{ "ADVECTION_DIFFUSION_2D" , PHYSICS::ADVECTION_DIFFUSION_2D , { "advection_diffusion_2D" }, "2D advection diffusion"},
		{ "ADVECTION_DIFFUSION_3D" , PHYSICS::ADVECTION_DIFFUSION_3D , { "advection_diffusion_3D" }, "3D advection diffusion"},
		{ "STRUCTURAL_MECHANICS_2D", PHYSICS::STRUCTURAL_MECHANICS_2D, { "structural_mechanics_2D" }, "2D structural mechanics"},
		{ "STRUCTURAL_MECHANICS_3D", PHYSICS::STRUCTURAL_MECHANICS_3D, { "structural_mechanics_3D" }, "3D structural mechanics"},
		{ "SHALLOW_WATER_2D"       , PHYSICS::SHALLOW_WATER_2D       , { "shallow_water_2D" }, "2D shallow water"}
	}));
};

}


#endif /* SRC_CONFIGURATION_GLOBALCONFIGURATION_H_ */
