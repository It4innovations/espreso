
#ifndef SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_
#define SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_

#include "config/description.h"

namespace espreso
{
struct PSOConfiguration : public ECFDescription
{
	int generations;
	double C1;
	double C2;
	double W_START;
	double W_END;

	PSOConfiguration();
};

struct DEConfiguration : public ECFDescription
{
	double F;
	double CR;

	DEConfiguration();
};

struct MicroDERConfiguration : public ECFDescription
{
	double F_START;
	double F_END;
	double CR;
	int RESTART_LIMIT;

	MicroDERConfiguration();
};

struct SOMAT3AConfiguration : public ECFDescription
{
	// TODO
};

struct SOMAConfiguration : public ECFDescription
{
	double PRT;
	double STEP;
	double PATH_LENGTH;

	SOMAConfiguration();
};

struct AutoOptimizationConfiguration : public ECFDescription
{
	enum class ALGORITHM
	{
		NONE,
		PARTICLE_SWARM,
		DIFFERENTIAL_EVOLUTION,
		SOMAT3A,
		SOMA,
		RANDOM,
		ALL_PERMUTATIONS,
		MICRO_DER
	};

	ALGORITHM algorithm;
	bool rounding_immediate;
	int population;

	PSOConfiguration particle_swarm;
	DEConfiguration differential_evolution;
	SOMAConfiguration soma;
	SOMAT3AConfiguration somat3a;
	MicroDERConfiguration micro_der;

	AutoOptimizationConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_ */
