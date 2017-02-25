
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
#include "physics/linearelasticity2d.h"
#include "physics/linearelasticity3d.h"
#include "physics/shallowwater2d.h"

namespace espreso {

enum class INPUT {
	WORKBENCH = 0,
	OPENFOAM = 1,
	ESDATA = 2,
	GENERATOR = 3
};

struct GlobalConfiguration: public Configuration {

	GlobalConfiguration(const std::string &file) { Reader::read(*this, file); Reader::set(this->env); }
	GlobalConfiguration(int *argc, char ***argv) { Reader::read(*this, argc, argv); Reader::set(this->env); }

	void print() { Reader::print(*this); }
	void store() { Reader::store(*this, { ".*" }); }

	OPTION(INPUT, input, "test input", INPUT::GENERATOR, OPTIONS({
			{ "WORKBENCH", INPUT::WORKBENCH, "Ansys Workbench input file" },
			{ "OPENFOAM", INPUT::OPENFOAM, "OpenFOAM input format" },
			{ "ESDATA", INPUT::ESDATA, "ESPRESO binary format" },
			{ "GENERATOR", INPUT::GENERATOR, "ESPRESO internal generator" }
	}));

	OPTION(PHYSICS, physics, "Used physics", PHYSICS::LINEAR_ELASTICITY_3D, OPTIONS({
		{ "LINEAR_ELASTICITY_2D"   , PHYSICS::LINEAR_ELASTICITY_2D   , "2D linear elasticity." },
		{ "LINEAR_ELASTICITY_3D"   , PHYSICS::LINEAR_ELASTICITY_3D   , "3D linear elasticity." },
		{ "TRANSIENT_ELASTICITY_2D", PHYSICS::TRANSIENT_ELASTICITY_2D, "2D transient elasticity." },
		{ "TRANSIENT_ELASTICITY_3D", PHYSICS::TRANSIENT_ELASTICITY_3D, "3D transient elasticity." },
		{ "ADVECTION_DIFFUSION_2D" , PHYSICS::ADVECTION_DIFFUSION_2D , "2D advection diffusion"},
		{ "ADVECTION_DIFFUSION_3D" , PHYSICS::ADVECTION_DIFFUSION_3D , "3D advection diffusion"},
		{ "SHALLOW_WATER_2D"       , PHYSICS::SHALLOW_WATER_2D       , "2D shallow water"},
		{ "STOKES"                 , PHYSICS::STOKES                 , "Stokes"}
	}));

	SUBCONFIG(Environment        , env         , "Environment dependent variables (set by ./env/threading.* scripts).");
	SUBCONFIG(OutputConfiguration, output      , "Output settings.");

	SUBCONFIG(ESPRESOGenerator   , generator   , "ESPRESO internal mesh generator.");
	SUBCONFIG(ESPRESOInput       , workbench   , "Mesh description in Ansys Workbench format.");
	SUBCONFIG(ESPRESOInput       , openfoam    , "Mesh description in OpenFOAM format.");
	SUBCONFIG(ESPRESOInput       , esdata      , "Mesh description in ESPRESO internal binary format.");

	SUBCONFIG(LinearElasticity2DConfiguration  , linear_elasticity_2D  , "2D Linear elasticity solver.");
	SUBCONFIG(LinearElasticity3DConfiguration  , linear_elasticity_3D  , "3D Linear elasticity solver.");
	SUBCONFIG(AdvectionDiffusion2DConfiguration, advection_diffusion_2D, "2D advection diffusiuon solver.");
	SUBCONFIG(AdvectionDiffusion3DConfiguration, advection_diffusion_3D, "3D advection diffusiuon solver.");
	SUBCONFIG(ShallowWater2DConfiguration      , shallow_water_2D      , "2D shallow water solver.");

	SUBCONFIG(Results, results, "Expected output results.");

	SUBCONFIG(Decomposer, decomposer, "./decomposer configuration");
};

}


#endif /* SRC_CONFIGURATION_GLOBALCONFIGURATION_H_ */
