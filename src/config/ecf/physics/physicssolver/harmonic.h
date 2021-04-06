
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_HARMONIC_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_HARMONIC_H_

#include "damping.h"
#include "config/description.h"

namespace espreso {

struct AlternatingFrequencyTime : public ECFDescription
{
	enum class TYPE {
		USER,
		AUTOMATIC
	};

	TYPE type;
	int time_samples;

	AlternatingFrequencyTime();
};

struct HarmonicSolverConfiguration: public ECFDescription {

	enum class INTERVAL_TYPE {
		LINEAR,
		LOGARITHMIC,
		OCTAVE_BAND,
//		USER
	};

	enum class OCTAVE_BAND {
		BAND_1,
		BAND_2,
		BAND_3,
		BAND_6,
		BAND_12,
		BAND_24,
	};

	INTERVAL_TYPE frequency_interval_type;

	double central_frequency;
	OCTAVE_BAND octave_band;

	double min_frequency, max_frequency;
	int num_samples;

	HarmonicDampingConfiguration damping;

	AlternatingFrequencyTime aft;

	HarmonicSolverConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_HARMONIC_H_ */
