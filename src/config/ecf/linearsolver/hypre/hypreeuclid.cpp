
#include "../../../configuration.hpp"
#include "hypreeuclid.h"

using namespace espreso;

HYPREEuclidConfiguration::HYPREEuclidConfiguration()
{
	levels = 1;
	REGISTER(levels, ECFMetaData()
			.setdescription({ "Set level k for ILU(k) factorization" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	set_bj = 0;
	REGISTER(set_bj, ECFMetaData()
			.setdescription({ "Use block Jacobi ILU preconditioning instead of PILU" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	stats = 0;
	REGISTER(stats, ECFMetaData()
			.setdescription({ "s summary of runtime settings and timing information" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	memory_stats = 0;
	REGISTER(memory_stats, ECFMetaData()
			.setdescription({ "a summary of Euclid’s memory usage" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	sparse_tol = 0.0;
	REGISTER(sparse_tol, ECFMetaData()
			.setdescription({ "Defines a drop tolerance for ILU(k)" })
			.setdatatype({ ECFDataType::FLOAT }));

	row_scale = 0;
	REGISTER(row_scale, ECFMetaData()
			.setdescription({ "values are scaled prior to factorization so that largest value in any row is +1 or -1" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));


	ilut_tol = 0.0;
	REGISTER(ilut_tol, ECFMetaData()
			.setdescription({ "uses ILUT and defines a drop tolerance relative to the largest absolute value of any entry in the row being factored)" })
			.setdatatype({ ECFDataType::FLOAT }));

}