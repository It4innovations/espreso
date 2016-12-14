
#ifndef SRC_CONFIG_GLOBALCONFIGURATION_H_
#define SRC_CONFIG_GLOBALCONFIGURATION_H_

#include "environment.h"
#include "output.h"

#include "input.h"
#include "inputgenerator.h"

#include "material.h"

#include "linearelasticity2d.h"
#include "linearelasticity3d.h"
#include "advectiondiffusion2d.h"
#include "advectiondiffusion3d.h"

#include "results.h"

#include "decomposer.h"

namespace espreso {

enum class INPUT {
	WORKBENCH = 0,
	OPENFOAM = 1,
	ESDATA = 2,
	GENERATOR = 3
};

enum class PHYSICS {
	LINEAR_ELASTICITY_2D,
	LINEAR_ELASTICITY_3D,
	TRANSIENT_ELASTICITY_2D,
	TRANSIENT_ELASTICITY_3D,
	ADVECTION_DIFFUSION_2D,
	ADVECTION_DIFFUSION_3D,
	STOKES
};

struct GlobalConfiguration: public Configuration {

	GlobalConfiguration(const std::string &file) { Reader::read(*this, file); Reader::set(*this); }
	GlobalConfiguration(int *argc, char ***argv) { Reader::read(*this, argc, argv); Reader::set(*this); }

	void print() { Reader::print(*this); }
	void store() { Reader::store(*this); }

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
		{ "STOKES"                 , PHYSICS::STOKES                 , "Stokes"}
	}));

	SUBCONFIG(Environment        , env         , "Environment dependent variables (set by ./env/threading.* scripts).");
	SUBCONFIG(OutputConfiguration, output      , "Output settings.");

	SUBCONFIG(ESPRESOGenerator   , generator   , "ESPRESO internal mesh generator.");
	SUBCONFIG(ESPRESOInput       , workbench   , "Mesh description in Ansys Workbench format.");
	SUBCONFIG(ESPRESOInput       , openfoam    , "Mesh description in OpenFOAM format.");
	SUBCONFIG(ESPRESOInput       , esdata      , "Mesh description in ESPRESO internal binary format.");
	SUBCONFIG(ESPRESOInput       , api         , "API description.");

	SUBCONFIG(LinearElasticity2DConfiguration  , linear_elasticity_2D  , "2D Linear elasticity solver.");
	SUBCONFIG(LinearElasticity3DConfiguration  , linear_elasticity_3D  , "3D Linear elasticity solver.");
	SUBCONFIG(AdvectionDiffusion2DConfiguration, advection_diffusion_2D, "2D advection diffusiuon solver.");
	SUBCONFIG(AdvectionDiffusion3DConfiguration, advection_diffusion_3D, "3D advection diffusiuon solver.");

	SUBCONFIG(Results, results, "Expected output results.");

	SUBCONFIG(Decomposer, decomposer, "./decomposer configuration");
};

}


#endif /* SRC_CONFIG_GLOBALCONFIGURATION_H_ */
