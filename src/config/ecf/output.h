
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

#include "config/description.h"
#include "physics/physics.h"

#include "physics/heattransfer.h"
#include "physics/structuralmechanics.h"

namespace espreso {

struct ECF;

struct MonitorConfiguration: public ECFDescription {

	enum class STATISTICS {
		MIN,
		MAX,
		AVG,
		NORM,
		ABSMIN,
		ABSMAX
	};

	std::string region;
	STATISTICS statistics;
	std::string property;

	MonitorConfiguration(const ECF *root);
protected:
	const ECF *_root;
};

struct ResultsSelectionConfiguration: public HeatTransferOutputSettings, public StructuralMechanicsOutputSettings {

	void basic()
	{
		HeatTransferOutputSettings::basic();
		StructuralMechanicsOutputSettings::basic();
		thickness = false;
	}

	void all()
	{
		HeatTransferOutputSettings::all();
		StructuralMechanicsOutputSettings::all();
		thickness = true;
	}

	bool thickness;

	ResultsSelectionConfiguration(const PhysicsConfiguration::TYPE &physics);
protected:
	const PhysicsConfiguration::TYPE &_physics;
};

struct HarmonicOuputConfiguration: public ECFDescription {

	enum class STORE_FREQUENCY {
		NEVER,
		EVERY_FREQUENCY,
		EVERY_NTH_FREQUENCY,
		SPECIFIC_FREQUENCIES
	};

	int samples;
	STORE_FREQUENCY results_store_frequency;

	std::vector<double> requested_frequencies; // TODO: implement array
	int results_nth_stepping;

	HarmonicOuputConfiguration();
};

struct OutputConfiguration: public ECFDescription {

	enum class FORMAT {
		VTK_LEGACY = 0,
		ENSIGHT,
		XDMF,
		STL_SURFACE,
		NETGEN
	};

	enum class WRITER {
		MPI,
		MPI_COLLECTIVE,
		POSIX
	};

	enum class LOGGER {
		USER,
		PARSER
	};

	enum class MODE {
		SYNC,
		PTHREAD,
	};

	enum class STORE_FREQUENCY {
		NEVER,
		EVERY_SUBSTEP,
		EVERY_NTH_SUBSTEP,
		LAST_SUBSTEP
	};

	enum class STORE_RESULTS {
		BASIC,
		ALL,
		USER
	};

	std::string log_dir;

	size_t verbose_level;
	size_t measure_level;
	std::string papi_events, papi_codes;
	LOGGER logger;

	size_t print_matrices, store_decomposition;

	FORMAT format;
	MODE mode;

	WRITER writer;
	size_t stripe_size, stripe_count;

//	size_t output_node_group_size;

	std::string path;

	STORE_FREQUENCY results_store_frequency, monitors_store_frequency;
	size_t results_nth_stepping, monitors_nth_stepping;

	STORE_RESULTS store_results;
	ResultsSelectionConfiguration results_selection;

	HarmonicOuputConfiguration frequency_to_time;

	bool settings, debug, catalyst;
	size_t catalyst_sleep_time;

	bool collected;

	std::map<size_t, MonitorConfiguration> monitoring;

	OutputConfiguration(const ECF *root);

protected:
	const ECF *_root;
};

}



#endif /* SRC_CONFIG_ECF_OUTPUT_H_ */
