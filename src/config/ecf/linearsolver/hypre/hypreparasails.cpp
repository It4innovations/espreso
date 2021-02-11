
#include "hypreparasails.h"

#include "config/configuration.hpp"

using namespace espreso;

HYPREParaSailsConfiguration::HYPREParaSailsConfiguration()
{
	threshold = 0.1;
	REGISTER(threshold, ECFMetaData()
			.setdescription({ "Set the threshold parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::FLOAT }));

	n_levels = 1;
	REGISTER(n_levels, ECFMetaData()
			.setdescription({ "Set the levels parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	filter = 0.1;
	REGISTER(filter, ECFMetaData()
			.setdescription({ "Set the filter parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::FLOAT }));

	symmetry = SYMMETRY::SPD;
	REGISTER(symmetry, ECFMetaData()
			.setdescription({ "Set the symmetry parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NON_INF").setdescription("nonsymmetric and/or indefinite problem, and nonsymmetric preconditioner"))
			.addoption(ECFOption().setname("SPD").setdescription("SPD problem, and SPD (factored) preconditioner"))
			.addoption(ECFOption().setname("NON_DEF_SPD").setdescription("nonsymmetric, definite problem, and SPD (factored) preconditioner")));

	loadbal = 0.0;
	REGISTER(loadbal, ECFMetaData()
			.setdescription({ "Set the load balance parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::FLOAT }));

	reuse = 0;
	REGISTER(reuse, ECFMetaData()
			.setdescription({ "Set the pattern reuse parameter for the ParaSails preconditioner" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));   

	logging = 3;
	REGISTER(logging, ECFMetaData()
			.setdescription({ "Set the logging parameter for the ParaSails preconditioner." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

}
