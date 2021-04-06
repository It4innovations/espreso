
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
		.setdescription({ "BoomerAMG settings." }));

	REGISTER(pcg, ECFMetaData()
		.setdescription({ "PCG settings." }));

	REGISTER(gmres, ECFMetaData()
		.setdescription({ "GMRES settings." }));

	REGISTER(flexgmres, ECFMetaData()
		.setdescription({ "FlexGMRES settings." }));

	REGISTER(lgmres, ECFMetaData()
		.setdescription({ "LGMRES settings." }));

	REGISTER(bicgstab, ECFMetaData()
		.setdescription({ "BiCGSTAB settings." }));

	REGISTER(cgnr, ECFMetaData()
		.setdescription({ "CGNR settings." }));

}



