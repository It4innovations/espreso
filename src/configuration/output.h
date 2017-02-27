
#ifndef SRC_CONFIGURATION_OUTPUT_H_
#define SRC_CONFIGURATION_OUTPUT_H_

#include "../configuration/configuration.hpp"

namespace espreso {

enum class OUTPUT_FORMAT {
	VTK_LEGACY_FORMAT = 0,
	VTK_BINARY_FORMAT = 1,
	VTK_MULTIBLOCK_FORMAT = 2,
	ENSIGHT_FORMAT = 3
};

struct OutputConfiguration: public Configuration {

	OPTION(OUTPUT_FORMAT, format, "Format - only LEGACY format is supported without VTK library", OUTPUT_FORMAT::VTK_LEGACY_FORMAT, OPTIONS({
		{ "VTK_LEGACY"    , OUTPUT_FORMAT::VTK_LEGACY_FORMAT    , "*.vtk files" },
		{ "VTK_BINARY"    , OUTPUT_FORMAT::VTK_BINARY_FORMAT    , "*.vtu files" },
		{ "VTK_MULTIBLOCK", OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT, "*.vtu + *.vtm files" },
		{ "ENSIGHT"       , OUTPUT_FORMAT::ENSIGHT_FORMAT       , "EnSight files" }
	}));

	PARAMETER(bool, compression, "Compression - needs VTK library", true);
	PARAMETER(double, decimation, "Decimation - needs VTK library", 0);

	PARAMETER(bool, results, "Save results", true);
	PARAMETER(bool, properties, "Save also input parameters", false);
	PARAMETER(bool, substeps, "Save substep results", false);
	PARAMETER(bool, gluing, "Save lagrange multipliers", false);

	PARAMETER(bool, catalyst, "Allow live visualization", false);
	PARAMETER(size_t, sleep, "Sleep interval between consecutive time step visialization", 0);

	PARAMETER(double, domain_shrink_ratio, "All domains are shrunk by this ratio", .95);
	PARAMETER(double, cluster_shrink_ratio  , "All clusters are shrunk by this ratio"  , .9);
};

}



#endif /* SRC_CONFIGURATION_OUTPUT_H_ */
