
#include "hypre.h"

#include "config/configuration.hpp"

using namespace espreso;

HYPREConfiguration::HYPREConfiguration()
{
	solver_type = SOLVER_TYPE::BoomerAMG;
	REGISTER(solver_type, ECFMetaData()
			.setdescription({ "Solver type...." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BoomerAMG").setdescription("Algebraic Multigrid"))
			.addoption(ECFOption().setname("PCG").setdescription("Preconditioned Conjugate Gradient Solver"))
			.addoption(ECFOption().setname("GMRES").setdescription("Generalized Minimal Residual Solver"))
			.addoption(ECFOption().setname("FlexGMRES").setdescription("Flexible Generalized Minimal Residual Solver"))
			.addoption(ECFOption().setname("LGMRES").setdescription("Accelerated Generalized Minimal Residual Solver"))
			.addoption(ECFOption().setname("BiCGSTAB").setdescription("Biconjugate Gradient Stabilized Solver "))
			.addoption(ECFOption().setname("CGNR").setdescription("Conjugate Gradient Method on the Normal Equations")));

	REGISTER(boomeramg, ECFMetaData()
		.setdescription({ "BoomerAMG settings." })
		.setcollapsed());

	REGISTER(pcg, ECFMetaData()
		.setdescription({ "PCG settings." })
		.setcollapsed());

	REGISTER(gmres, ECFMetaData()
		.setdescription({ "GMRES settings." })
		.setcollapsed());

	REGISTER(flexgmres, ECFMetaData()
		.setdescription({ "FlexGMRES settings." })
		.setcollapsed());

	REGISTER(lgmres, ECFMetaData()
		.setdescription({ "LGMRES settings." })
		.setcollapsed());

	REGISTER(bicgstab, ECFMetaData()
		.setdescription({ "BiCGSTAB settings." })
		.setcollapsed());

	REGISTER(cgnr, ECFMetaData()
		.setdescription({ "CGNR settings." })
		.setcollapsed());

}



