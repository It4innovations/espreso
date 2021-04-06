
#include "hyprebicgstab.h"

#include "config/configuration.hpp"

using namespace espreso;

HYPREBiCGSTABConfiguration::HYPREBiCGSTABConfiguration()
{
	relative_conv_tol = 1e-8;
	REGISTER(relative_conv_tol, ECFMetaData()
			.setdescription({ "Set the relative convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

	absolute_conv_tol = 0;
	REGISTER(absolute_conv_tol, ECFMetaData()
			.setdescription({ "Set the absolute convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

	max_iterations = 1000;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Set maximum number of iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	preconditioner = PRECONDITIONER::BoomerAMG;
	REGISTER(preconditioner, ECFMetaData()
			.setdescription({ "Preconditioner" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BoomerAMG").setdescription("Set BoomerAMG as a preconditioner"))
			.addoption(ECFOption().setname("ParaSails").setdescription("Set ParaSails as a preconditioner"))
			.addoption(ECFOption().setname("Euclid").setdescription("Set Euclid as a preconditioner"))
			.addoption(ECFOption().setname("Pilut").setdescription("Set Pilut as a preconditioner"))
			.addoption(ECFOption().setname("NONE").setdescription("Solver without preconditioner")));

	REGISTER(boomeramg, ECFMetaData()
			.setdescription({ "BoomerAMG settings." }));

	REGISTER(parasails, ECFMetaData()
			.setdescription({ "ParaSails settings." }));

	REGISTER(euclid, ECFMetaData()
			.setdescription({ "Euclid settings." }));

	REGISTER(pilut, ECFMetaData()
			.setdescription({ "Pilut settings." }));

	solver_info = SOLVER_INFO::NO_INFO;
	REGISTER(solver_info, ECFMetaData()
			.setdescription({ "Print solver info" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NO_INFO").setdescription("no printout"))
			.addoption(ECFOption().setname("SETUP_INFO").setdescription("print setup information"))
			.addoption(ECFOption().setname("SOLVE_INFO").setdescription("print solve information"))
			.addoption(ECFOption().setname("SETUP_SOLVE_INFO").setdescription("print both setup and solve information")));

}
