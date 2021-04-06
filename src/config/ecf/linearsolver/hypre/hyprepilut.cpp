
#include "../../../configuration.hpp"
#include "hyprepilut.h"

using namespace espreso;

HYPREPilutConfiguration::HYPREPilutConfiguration()
{
	max_iter = 1;
	REGISTER(max_iter, ECFMetaData()
			.setdescription({ "Set maximum number of iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	drop_tol = 0.0001;
	REGISTER(drop_tol, ECFMetaData()
			.setdescription({ "Set maximum number of iterations" })
			.setdatatype({ ECFDataType::FLOAT }));

	row_size = 20;
	REGISTER(row_size, ECFMetaData()
			.setdescription({ "Set the maximum number of nonzeros to be retained in each row of L and U" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));


}

