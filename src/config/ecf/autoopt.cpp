
#include "autoopt.h"
#include "config/configuration.hpp"

using namespace espreso;

PSOConfiguration::PSOConfiguration()
{
	this->generations = 10;
	REGISTER(generations, ECFMetaData()
		.setdescription({"Number of migrations"})
		.setdatatype({ECFDataType::INTEGER}));

	this->C1 = 2.0f;
	REGISTER(C1, ECFMetaData()
		.setdescription({"C1"})
		.setdatatype({ECFDataType::FLOAT}));

	this->C2 = 2.0f;
	REGISTER(C2, ECFMetaData()
		.setdescription({"C2"})
		.setdatatype({ECFDataType::FLOAT}));

	this->W_START = 0.9f;
	REGISTER(W_START, ECFMetaData()
		.setdescription({"W_START"})
		.setdatatype({ECFDataType::FLOAT}));

	this->W_END = 0.4f;
	REGISTER(W_END, ECFMetaData()
		.setdescription({"W_END"})
		.setdatatype({ECFDataType::FLOAT}));
}

DEConfiguration::DEConfiguration()
{
	this->F = 0.9f;
	REGISTER(F, ECFMetaData()
		.setdescription({"F"})
		.setdatatype({ECFDataType::FLOAT}));

	this->CR = 0.4f;
	REGISTER(CR, ECFMetaData()
		.setdescription({"CR"})
		.setdatatype({ECFDataType::FLOAT}));
}

MicroDERConfiguration::MicroDERConfiguration()
{
	this->F_START = 0.1f;
	REGISTER(F_START, ECFMetaData()
		.setdescription({"F_START"})
		.setdatatype({ECFDataType::FLOAT}));

	this->F_END = 1.5f;
	REGISTER(F_END, ECFMetaData()
		.setdescription({"F_END"})
		.setdatatype({ECFDataType::FLOAT}));

	this->CR = 0.9f;
	REGISTER(CR, ECFMetaData()
		.setdescription({"CR"})
		.setdatatype({ECFDataType::FLOAT}));

	this->RESTART_LIMIT = 5;
	REGISTER(RESTART_LIMIT, ECFMetaData()
		.setdescription({"RESTART_LIMIT"})
		.setdatatype({ECFDataType::NONNEGATIVE_INTEGER}));
}

SOMAConfiguration::SOMAConfiguration()
{
	this->PRT = 0.3f;
	REGISTER(PRT, ECFMetaData()
		.setdescription({"PRT"})
		.setdatatype({ECFDataType::FLOAT}));

	this->STEP = 0.8f;
	REGISTER(STEP, ECFMetaData()
		.setdescription({"STEP"})
		.setdatatype({ECFDataType::FLOAT}));

	this->PATH_LENGTH = 3.0f;
	REGISTER(PATH_LENGTH, ECFMetaData()
		.setdescription({"PATH_LENGTH"})
		.setdatatype({ECFDataType::FLOAT}));
}

AutoOptimizationConfiguration::AutoOptimizationConfiguration()
{
	this->algorithm = ALGORITHM::NONE;
	REGISTER(algorithm, ECFMetaData()
		.setdescription({"Biologically inspired algorithm"})
		.setdatatype({ECFDataType::OPTION})
		.addoption(ECFOption().setname("NONE").setdescription("Automatic optimization is not used."))
		.addoption(ECFOption().setname("PARTICLE_SWARM").setdescription("Particle Swarm Optimization."))
		.addoption(ECFOption().setname("DIFFERENTIAL_EVOLUTION").setdescription("Differential Evolution."))
		.addoption(ECFOption().setname("SOMAT3A").setdescription("SOMAT3A."))
		.addoption(ECFOption().setname("SOMA").setdescription("SOMA."))
		.addoption(ECFOption().setname("RANDOM").setdescription("Generating random configurations."))
		.addoption(ECFOption().setname("ALL_PERMUTATIONS").setdescription("All parameter permutations."))
		.addoption(ECFOption().setname("MICRO_DER").setdescription("Micro Differential Evolution Ray.")));

	this->rounding_immediate = false;
	REGISTER(rounding_immediate, ECFMetaData()
		.setdescription({"Immediate rounding"})
		.setdatatype({ECFDataType::BOOL}));

	this->population = 10;
	REGISTER(population, ECFMetaData()
		.setdescription({"Population"})
		.setdatatype({ECFDataType::INTEGER}));

	REGISTER(particle_swarm, ECFMetaData()
		.setdescription({"Particle Swarm Optimization"}));
	REGISTER(differential_evolution, ECFMetaData()
		.setdescription({"Differential Evolution"}));
	REGISTER(soma, ECFMetaData()
		.setdescription({"SOMA"}));
	REGISTER(micro_der, ECFMetaData()
		.setdescription({"Micro Differential Evolution Ray"}));
	// REGISTER(somat3a, ECFMetaData()
	// 	.setdescription({"SOMA Team-to-Team"}));
}
