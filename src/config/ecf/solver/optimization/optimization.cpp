
#include "optimization.h"
#include "../../../configuration.hpp"

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

OptimizationConfiguration::OptimizationConfiguration()
{
	this->algorithm = ALGORITHM::PARTICLE_SWARM;
	REGISTER(algorithm, ECFMetaData()
		.setdescription({"Biologically inspired algorithm"})
		.setdatatype({ECFDataType::OPTION})
		.addoption(ECFOption().setname("PARTICLE_SWARM").setdescription("Particle Swarm Optimization"))
		.addoption(ECFOption().setname("DIFFERENTIAL_EVOLUTION").setdescription("Differential Evolution"))
		.addoption(ECFOption().setname("SOMAT3A").setdescription("SOMAT3A")));
	
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
	// REGISTER(somat3a, ECFMetaData()
	// 	.setdescription({"SOMA Team-to-Team"}));
}