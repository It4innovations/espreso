
#ifndef SRC_CONFIGURATION_OUTPUT_H_
#define SRC_CONFIGURATION_OUTPUT_H_

#include "../configuration/configuration.hpp"

namespace espreso {

enum class OUTPUT_FORMAT {
	VTK_LEGACY = 0,
	VTK_XML_ASCII = 1,
	VTK_XML_BINARY = 2,
	ENSIGHT = 3
};

enum class OUTPUT_MODE {
	SYNC,
	THREAD,
	MPI,
};

struct OutputConfiguration: public Configuration {

	OPTION(OUTPUT_FORMAT, format, "Format - only LEGACY format is supported without VTK library", OUTPUT_FORMAT::VTK_XML_ASCII, OPTIONS({
		{ "VTK_LEGACY"    , OUTPUT_FORMAT::VTK_LEGACY    , "*.vtk files" },
		{ "VTK_XML_ASCII" , OUTPUT_FORMAT::VTK_XML_ASCII , "*.vtu files in ASCII format" },
		{ "VTK_XML_BINARY", OUTPUT_FORMAT::VTK_XML_BINARY, "*.vtu files in binary format" },
		{ "ENSIGHT"       , OUTPUT_FORMAT::ENSIGHT       , "EnSight files" }
	}));

	OPTION(OUTPUT_MODE, mode, "Mode of ASYNC library", OUTPUT_MODE::THREAD, OPTIONS({
		{ "SYNC"  , OUTPUT_MODE::SYNC  , "Storing is synchronized." },
		{ "THREAD", OUTPUT_MODE::THREAD, "Storing is performed by the last thread." },
		{ "MPI"   , OUTPUT_MODE::MPI   , "Storing is forwarded to I/O MPI nodes according to OUTPUT_NODE_GROUP_SIZE." },
	}));

	PARAMETER(size_t, output_node_group_size, "Max number of compution nodes for one output node.", 7);

	PARAMETER(std::string, path, "Path to output files.", "results");

	PARAMETER(bool, compression, "Compression - needs VTK library", true);
	PARAMETER(double, decimation, "Decimation - needs VTK library", 0);

	PARAMETER(bool, results, "Save results", true);
	PARAMETER(bool, settings, "Save also input parameters", false);
	PARAMETER(bool, iterations, "Save results for all iterations", false);
	PARAMETER(bool, FETI_data, "Save FETI data (decomposition, fix points, gluing, etc...)", false);

	PARAMETER(bool, catalyst, "Allow live visualization", false);
	PARAMETER(size_t, sleep, "Sleep interval between consecutive time step visialization", 0);

	PARAMETER(bool  , collected             , "Gather results from all processes to one file." , false);
	PARAMETER(bool  , separate_bodies       , "Store mesh bodies to separated files."          , true);
	PARAMETER(bool  , separate_materials    , "Store mesh material to separated files."        , false);

	PARAMETER(double, domain_shrink_ratio   , "All domains are shrunk by this ratio (effective only for COLLECTED=FALSE)."   , .95);
	PARAMETER(double, cluster_shrink_ratio  , "All clusters are shrunk by this ratio (effective only for COLLECTED=FALSE)."  , .9);

	SUBMULTIMAP(std::string, std::string, monitoring, "Results statistics in some regions. OPERATION = { AVERAGE, MIN, MAX }", "REGION", "<OPERATION> <VARIABLE>");
};

}



#endif /* SRC_CONFIGURATION_OUTPUT_H_ */
