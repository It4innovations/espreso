
#include "input.h"

#include "config/configuration.hpp"

using namespace espreso;

InputConfiguration::InputConfiguration()
{
	path = ".";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Path" })
			.setdatatype({ ECFDataType::STRING }));

	format = FORMAT::ANSYS_CDB;
	REGISTER(format, ECFMetaData()
			.setdescription({ "An input data type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ANSYS_CDB").setdescription("Ansys WorkCDB format."))
			.addoption(ECFOption().setname("OpenFOAM").setdescription("OpenFOAM format."))
			.addoption(ECFOption().setname("ABAQUS").setdescription("ABAQUS format."))
			.addoption(ECFOption().setname("XDMF").setdescription("XDMF format."))
			.addoption(ECFOption().setname("ENSIGHT").setdescription("Ensigh Gold format."))
			.addoption(ECFOption().setname("VTK_LEGACY").setdescription("VTK Legacy format."))
			.addoption(ECFOption().setname("NETGEN").setdescription("Neutral Netgen format.")));

	keep_material_sets = false;
	REGISTER(keep_material_sets, ECFMetaData()
			.setdescription({ "Keep material sets" })
			.setdatatype({ ECFDataType::BOOL }));

	convert_database = false;
	REGISTER(convert_database, ECFMetaData()
			.setdescription({ "Convert database" })
			.setdatatype({ ECFDataType::BOOL }));

	duplication_tolerance = 1e-8;
	REGISTER(duplication_tolerance, ECFMetaData()
			.setdescription({ "Tolerance for merging nodes according to coordinates." })
			.setdatatype({ ECFDataType::FLOAT }));

	loader = LOADER::MPI;
	REGISTER(loader, ECFMetaData()
			.setdescription({ "A type of used function for loading data." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("MPI").setdescription("Use MPI_File_read_at."))
			.addoption(ECFOption().setname("MPI_COLLECTIVE").setdescription("Use MPI_File_read_at_all."))
			.addoption(ECFOption().setname("POSIX").setdescription("Use POSIX API.")));

	stripe_size = 1024 * 1024;
	REGISTER(stripe_size, ECFMetaData()
			.setdescription({ "Stripe size of the input file." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	third_party_scalability_limit = 768;
	REGISTER(third_party_scalability_limit, ECFMetaData()
			.setdescription({ "Maximum number of MPI processes used for non-scalable routines." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	ecfdescription->addSeparator();

	REGISTER(transformations, ECFMetaData()
			.setdescription({ "List of transformations", "Transformation" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "ALL_ELEMENTS" }));

	REGISTER(decomposition, ECFMetaData()
			.setdescription({ "Domains decomposition settings." }));
}

