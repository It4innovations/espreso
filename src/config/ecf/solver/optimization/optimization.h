
#ifndef SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_
#define SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_

#include "../../../configuration.h"

namespace espreso 
{
struct PSOConfiguration : public ECFObject
{
	int generations;
	double C1;
	double C2;
	double W_START;
	double W_END;

	PSOConfiguration();
};

struct DEConfiguration : public ECFObject
{
	double F;
	double CR;

	DEConfiguration();
};


struct SOMAT3AConfiguration : public ECFObject
{
	// TODO
};

struct OptimizationConfiguration : public ECFObject
{
	enum class ALGORITHM
	{
		PARTICLE_SWARM,
    	DIFFERENTIAL_EVOLUTION,
    	SOMAT3A
	};

	ALGORITHM algorithm;
	bool rounding_immediate;
	int population;

	PSOConfiguration particle_swarm;
	DEConfiguration differential_evolution;
	SOMAT3AConfiguration somat3a;

	OptimizationConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_SOLVER_OPTIMIZATION_H_ */