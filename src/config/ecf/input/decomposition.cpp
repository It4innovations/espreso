
#include "decomposition.h"
#include "config/configuration.hpp"

using namespace espreso;

ParMETISConfiguration::ParMETISConfiguration()
{
	refinement = false;
	REGISTER(refinement, ECFMetaData()
			.setdescription({ "Apply refinement." })
			.setdatatype({ ECFDataType::BOOL }));

	tolerance = 1.05;
	REGISTER(tolerance, ECFMetaData()
			.setdescription({ "Imbalance tolerance." })
			.setdatatype({ ECFDataType::FLOAT }));
}

METISConfiguration::METISConfiguration()
{
	objective_type = OBJECTIVE_TYPE::VOLUME;
	REGISTER(objective_type, ECFMetaData()
			.setdescription({ "Decomposition respect materials." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VOLUME").setdescription("METIS tries to minimize communication volume. Recommended to examples with various element size."))
			.addoption(ECFOption().setname("EDGECUT").setdescription("METIS tries to minimize edgecut.")));

	continuous = 1;
}

PTScotchConfiguration::PTScotchConfiguration()
{

}

ScotchConfiguration::ScotchConfiguration()
{

}

KaHIPConfiguration::KaHIPConfiguration() {

}

DecompositionConfiguration::DecompositionConfiguration()
{
	parallel_decomposer = ParallelDecomposer::PARMETIS;
	REGISTER(parallel_decomposer, ECFMetaData()
			.setdescription({ "Tool that is used for decomoposition. " })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Skip parallel decomposition."))
			.addoption(ECFOption().setname("METIS").setdescription("METIS library (called by the root process)."))
			.addoption(ECFOption().setname("PARMETIS").setdescription("ParMETIS library."))
			.addoption(ECFOption().setname("PTSCOTCH").setdescription("PT-Scotch library."))
			.addoption(ECFOption().setname("HILBERT_CURVE").setdescription("Sort elements centers according Hilbert space filling curve.")));

	sequential_decomposer = SequentialDecomposer::METIS;
	REGISTER(sequential_decomposer, ECFMetaData()
			.setdescription({ "Tool that is used for decomoposition. " })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Skip sequential decomposition."))
			.addoption(ECFOption().setname("METIS").setdescription("METIS library."))
			.addoption(ECFOption().setname("Scotch").setdescription("Scotch library."))
			.addoption(ECFOption().setname("KaHIP").setdescription("KaHIP library.")));

	mesh_duplication = 1;
	REGISTER(mesh_duplication, ECFMetaData()
			.setdescription({ "The number of parallel mesh instances (usefull etc. for harmonic solvers)." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	domains = 0;
	REGISTER(domains, ECFMetaData()
			.setdescription({ "Number of domains for each cluster (Keep 0 for automatic decomposition)." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.setautooptrange(1, 512));

	ecfdescription->addSpace();

	force_continuity = false;
	REGISTER(force_continuity, ECFMetaData()
			.setdescription({ "Force continuous decomposition." })
			.setdatatype({ ECFDataType::BOOL }));

	separate_materials = false;
	separate_regions = false;
	separate_etypes = false;
	REGISTER(separate_materials, ECFMetaData()
			.setdescription({ "Decomposition respect materials." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_regions, ECFMetaData()
			.setdescription({ "Decomposition respect regions." })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(separate_etypes, ECFMetaData()
			.setdescription({ "Decomposition respect elements types." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(parmetis_options, ECFMetaData()
			.setdescription({ "ParMETIS options." }));
	REGISTER(metis_options, ECFMetaData()
			.setdescription({ "ParMETIS options." }));
	REGISTER(ptscotch_options, ECFMetaData()
			.setdescription({ "PT-Scotch options." }));
	REGISTER(scotch_options, ECFMetaData()
			.setdescription({ "PT-Scotch options." }));
	REGISTER(kahip_options, ECFMetaData()
			.setdescription({ "KaHIP options." }));
}




