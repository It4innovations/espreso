
#ifndef SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_
#define SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct ParMETISConfiguration: public ECFDescription {

	bool refinement;

	ParMETISConfiguration();
};

struct METISConfiguration: public ECFDescription {

	enum class OBJECTIVE_TYPE {
		VOLUME, EDGECUT
	};

	OBJECTIVE_TYPE objective_type;
	int continuous;

	METISConfiguration();
};

struct PTScotchConfiguration: public ECFDescription {

	PTScotchConfiguration();
};

struct ScotchConfiguration: public ECFDescription {

	ScotchConfiguration();
};

struct KaHIPConfiguration: public ECFDescription {

	KaHIPConfiguration();
};

struct DecompositionConfiguration: public ECFDescription {

	enum class ParallelDecomposer {
		NONE,
		METIS,
		PARMETIS,
		PTSCOTCH,
		HILBERT_CURVE
	};

	enum class SequentialDecomposer {
		NONE,
		METIS,
		SCOTCH,
		KAHIP
	};

	ParallelDecomposer parallel_decomposer;
	SequentialDecomposer sequential_decomposer;
	int mesh_duplication;
	int domains;

	bool separate_materials, separate_regions, separate_etypes;
	ParMETISConfiguration parmetis_options;
	METISConfiguration metis_options;
	PTScotchConfiguration ptscotch_options;
	ScotchConfiguration scotch_options;
	KaHIPConfiguration kahip_options;

	DecompositionConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_DECOMPOSITION_H_ */
