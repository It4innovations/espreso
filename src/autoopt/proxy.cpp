
#include "proxy.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/output.h"
#include "config/ecf/autoopt.h"
#include "wrappers/mpi/communication.h"

#include <limits>
#include <chrono>
#include <sstream>

using namespace espreso;

ParameterManager::ParameterManager(std::vector<ECFParameter*>& parameters, int population, bool roundingImmediate)
: m_params(parameters), m_immediate(roundingImmediate),
  m_generator(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
  m_dist_decimal(0.0f, 1.0f), m_dist_dimension(0, parameters.size() - 1),
  m_dist_population(0, population - 1)
{
	this->m_verconfs = {
		{0, 2, 0, 0, 0},
		{0, 2, 0, 1, 0},
		{0, 2, 1, 0, 0},
		{0, 2, 1, 1, 0},
		{0, 3, 0, 0, 0},
		{0, 3, 0, 1, 0},
		{0, 3, 1, 0, 0},
		{0, 3, 1, 1, 0},
		{1, 0, 1, 0, 0},
		{1, 0, 1, 1, 0},
		{1, 1, 0, 0, 0},
		{1, 2, 0, 0, 0},
		{1, 2, 0, 1, 0},
		{1, 2, 1, 0, 0},
		{1, 2, 1, 1, 0},
		{1, 3, 0, 0, 0},
		{1, 3, 0, 1, 0},
		{1, 3, 1, 0, 0},
		{1, 3, 1, 1, 0},
		{1, 4, 0, 0, 0},
		{1, 4, 0, 1, 0},
		{1, 4, 1, 0, 0},
		{1, 4, 1, 1, 0},
		{2, 2, 0, 0, 0},
		{2, 2, 0, 1, 0},
		{2, 2, 1, 0, 0},
		{2, 2, 1, 1, 0},
		{2, 3, 0, 0, 0},
		{2, 3, 0, 1, 0},
		{2, 3, 1, 0, 0},
		{2, 3, 1, 1, 0},
		{3, 0, 0, 0, 0},
		{3, 0, 0, 1, 0},
		{3, 0, 1, 1, 0},
		{3, 1, 0, 0, 0},
		{3, 1, 1, 0, 0},
		{3, 1, 1, 1, 0},
		{3, 2, 0, 0, 0},
		{3, 2, 0, 1, 0},
		{3, 2, 1, 0, 0},
		{3, 2, 1, 1, 0},
		{3, 3, 0, 0, 0},
		{3, 3, 0, 1, 0},
		{3, 3, 1, 0, 0},
		{3, 3, 1, 1, 0},
		{3, 4, 0, 0, 0},
		{3, 4, 0, 1, 0},
		{3, 4, 1, 0, 0},
		{3, 4, 1, 1, 0},
		{3, 6, 0, 0, 0},
		{3, 6, 0, 1, 0},
		{3, 6, 1, 0, 0},
		{3, 6, 1, 1, 0},
	};
}

int ParameterManager::count() const
{
	return m_params.size();
}

int ParameterManager::generateInt(int start, int endInclusive)
{
	std::uniform_int_distribution<int> d(start, endInclusive);
	return d(m_generator);
}

double ParameterManager::generateDouble(double start, double endInclusive)
{
	std::uniform_real_distribution<double> d(start, endInclusive);
	return d(m_generator);
}

int ParameterManager::generateSpecimenNumber()
{
	return m_dist_population(m_generator);
}

int ParameterManager::generateParameterNumber()
{
	return m_dist_dimension(m_generator);
}

double ParameterManager::generateDecimal()
{
	return m_dist_decimal(m_generator);
}
// double ParameterManager::generateDouble()
// {
//	 return static_cast<double>(m_generator()) / static_cast<double>(m_generator.max());
// }

std::vector<double> ParameterManager::generateConfiguration()
{
	return this->generateRandomConfiguration();
}

double ParameterManager::generateParameter(int id)
{
	double value;
	ECFDataType dt = m_params[id]->metadata.datatype.front();
	
	if (dt == ECFDataType::INTEGER)
	{
		value = m_generator() - m_generator() ;
	}
	else if (dt == ECFDataType::POSITIVE_INTEGER)
	{
		value = m_generator() % m_generator.max() + 1;
	}
	else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
	{
		if (m_params[id]->metadata.range)
		{
			std::uniform_int_distribution<int> dist(
				std::atoi(m_params[id]->metadata.range->min.c_str()),
				std::atoi(m_params[id]->metadata.range->max.c_str())
			);
			value = dist(m_generator);
		}
		else 
		{ value = m_generator(); }
	}
	else if (dt == ECFDataType::FLOAT)
	{
		value = (double)(m_generator() - m_generator()) + generateDecimal();
	}
	else if (dt == ECFDataType::OPTION)
	{
		std::uniform_int_distribution<int> dist(0, m_params[id]->metadata.options.size() - 1);
		value = dist(m_generator);
	}
	else if (dt == ECFDataType::ENUM_FLAGS)
	{
		std::uniform_int_distribution<int> dist(0, m_params[id]->metadata.options.size() - 1);
		value = dist(m_generator);
	}
	else if (dt == ECFDataType::BOOL)
	{
		std::uniform_int_distribution<int> dist(0, 1);
		value = dist(m_generator);
	}
	else
	{
		eslog::globalerror("Optimization - UNKNOWN PARAMETER DATA TYPE!\n");
	}

	return value;
}

std::vector<double> ParameterManager::generateRandomConfiguration()
{
	std::vector<double> configuration;
	for (auto it = m_params.begin(); it != m_params.end(); ++it)
	{
		ECFDataType dt = (*it)->metadata.datatype.front();
		if (dt == ECFDataType::INTEGER)
		{
			configuration.push_back( m_generator() - m_generator() );
		}
		else if (dt == ECFDataType::POSITIVE_INTEGER)
		{
			configuration.push_back(m_generator() % m_generator.max() + 1);
		}
		else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
		{
			if ((*it)->metadata.range)
			{
				std::uniform_int_distribution<int> dist(
					std::atoi((*it)->metadata.range->min.c_str()),
					std::atoi((*it)->metadata.range->max.c_str())
				);
				configuration.push_back( dist(m_generator) );
			}
			else 
			{ configuration.push_back(m_generator()); }
		}
		else if (dt == ECFDataType::FLOAT)
		{
			configuration.push_back(
				(double)(m_generator() - m_generator())
				+ generateDecimal()
			);
		}
		else if (dt == ECFDataType::OPTION)
		{
			std::uniform_int_distribution<int> dist(0, (*it)->metadata.options.size() - 1);
			configuration.push_back( dist(m_generator) );
		}
		else if (dt == ECFDataType::ENUM_FLAGS)
		{
			std::uniform_int_distribution<int> dist(0, (*it)->metadata.options.size() - 1);
			configuration.push_back( dist(m_generator) );
		}
		else if (dt == ECFDataType::BOOL)
		{
			std::uniform_int_distribution<int> dist(0, 1);
			configuration.push_back( dist(m_generator) );
		}
		else
		{
			eslog::globalerror("Optimization - UNKNOWN PARAMETER DATA TYPE!\n");
		}
	}

	return configuration;
}

std::vector<double> ParameterManager::generateVerifiedConfiguration()
{
	std::vector<double> configuration;
	
	std::uniform_int_distribution<int> dist(0, m_verconfs.size() - 1);
	std::vector<double>& verconf = this->m_verconfs[dist(m_generator)];

	for (auto it = m_params.begin(); it != m_params.end(); ++it)
	{
		ECFDataType dt = (*it)->metadata.datatype.front();
		
		if ((*it)->name.compare("preconditioner") == 0)
		{ configuration.push_back(verconf[0]); }
		else if ((*it)->name.compare("iterative_solver") == 0)
		{ configuration.push_back(verconf[1]); }
		else if ((*it)->name.compare("redundant_lagrange") == 0)
		{ configuration.push_back(verconf[2]); }
		else if ((*it)->name.compare("scaling") == 0)
		{ configuration.push_back(verconf[3]); }
		else if ((*it)->name.compare("method") == 0)
		{ configuration.push_back(verconf[4]); }
		else if (dt == ECFDataType::INTEGER)
		{
			configuration.push_back( m_generator() - m_generator() );
		}
		else if (dt == ECFDataType::POSITIVE_INTEGER)
		{
			configuration.push_back(m_generator() % m_generator.max() + 1);
		}
		else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
		{
			if ((*it)->metadata.range)
			{
				std::uniform_int_distribution<int> dist(
					std::atoi((*it)->metadata.range->min.c_str()),
					std::atoi((*it)->metadata.range->max.c_str())
				);
				configuration.push_back( dist(m_generator) );
			}
			else 
			{ configuration.push_back(m_generator()); }
		}
		else if (dt == ECFDataType::FLOAT)
		{
			configuration.push_back(
				(double)(m_generator() - m_generator())
				+ generateDecimal()
			);
		}
		else if (dt == ECFDataType::OPTION)
		{
			std::uniform_int_distribution<int> dist(0, (*it)->metadata.options.size() - 1);
			configuration.push_back( dist(m_generator) );
		}
		else if (dt == ECFDataType::ENUM_FLAGS)
		{
			std::uniform_int_distribution<int> dist(0, (*it)->metadata.options.size() - 1);
			configuration.push_back( dist(m_generator) );
		}
		else if (dt == ECFDataType::BOOL)
		{
			std::uniform_int_distribution<int> dist(0, 1);
			configuration.push_back( dist(m_generator) );
		}
		else
		{
			eslog::globalerror("Optimization - UNKNOWN PARAMETER DATA TYPE!\n");
		}
	}

	return configuration;
}

double ParameterManager::getParameterMin(int id)
{
	ECFParameter* p = this->m_params[id];
	ECFDataType dt = p->metadata.datatype.front();

	if (dt == ECFDataType::INTEGER)
	{
		return std::numeric_limits<int>::min();
	}
	else if (dt == ECFDataType::POSITIVE_INTEGER)
	{
		return 1;
	}
	else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
	{
		if (p->metadata.range)
		{ return std::atoi(p->metadata.range->min.c_str()); }
		else 
		{ return 0; }
	}
	else if (dt == ECFDataType::OPTION)
	{
		return 0;
	}
	else if (dt == ECFDataType::ENUM_FLAGS)
	{
		return 0;
	}
	else if (dt == ECFDataType::BOOL)
	{
		return 0;
	}
	else
	{
		eslog::globalerror("Optimizer: Unknown parameter\n");
		return 0;
	}
}

double ParameterManager::getParameterMax(int id)
{
	ECFParameter* p = this->m_params[id];
	ECFDataType dt = p->metadata.datatype.front();

	if (dt == ECFDataType::INTEGER)
	{
		return std::numeric_limits<int>::max();
	}
	else if (dt == ECFDataType::POSITIVE_INTEGER)
	{
		return std::numeric_limits<int>::max();
	}
	else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
	{
		if (p->metadata.range)
		{ return std::atoi(p->metadata.range->max.c_str()); }
		else 
		{ return std::numeric_limits<int>::max(); }
	}
	else if (dt == ECFDataType::OPTION)
	{
		int size = p->metadata.options.size();
		return size - 1;
	}
	else if (dt == ECFDataType::ENUM_FLAGS)
	{
		int size = p->metadata.options.size();
		return size - 1;
	}
	else if (dt == ECFDataType::BOOL)
	{
		return 1;
	}
	else
	{
		eslog::globalerror("Optimizer: Unknown parameter\n");
		return 0;
	}
}

double ParameterManager::checkParameter(int id, double value)
{
	if (m_immediate) return this->_checkParameter_immediate_rounding(id, value);
	else return this->_checkParameter_no_rounding(id, value);
}

bool ParameterManager::areParameterValuesEqual(int id, double val1, double val2)
{
	ECFParameter* p = this->m_params[id];
	ECFDataType dt = p->metadata.datatype.front();

	switch (dt) {
		case ECFDataType::INTEGER:
		case ECFDataType::POSITIVE_INTEGER:
		case ECFDataType::NONNEGATIVE_INTEGER:
		case ECFDataType::ENUM_FLAGS:
		case ECFDataType::OPTION:
		case ECFDataType::BOOL: {
			int ival1 = (int)val1;
			int ival2 = (int)val2;
			if (ival1 == ival2) return true;
			else return false;
		}
		default:
			eslog::globalerror("Optimizer: Comparing values of parameters with unknown datatype!");
			return false;
	}
}

double ParameterManager::_checkParameter_immediate_rounding(int id, double value)
{
	ECFParameter* p = this->m_params[id];
	ECFDataType dt = p->metadata.datatype.front();

	if (dt == ECFDataType::INTEGER)
	{
		return (int)value;
	}
	else if (dt == ECFDataType::POSITIVE_INTEGER)
	{
		int val = (int)value;
		if (val < 1) return 1;
		return val;
	}
	else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
	{
		int val = (int)value;
		if (p->metadata.range)
		{
			if (val < std::atoi(p->metadata.range->min.c_str()))
			{ return std::atoi(p->metadata.range->min.c_str()); }
			else if (val > std::atoi(p->metadata.range->max.c_str()))
			{ return std::atoi(p->metadata.range->max.c_str()); }
		}
		else if (val < 0) return 0;
	}
	else if (dt == ECFDataType::OPTION)
	{
		int opt = (int)value;
		if (opt < 0) return 0;
		int size = p->metadata.options.size();
		if (opt >= size) return size - 1;
		return opt;
	}
	else if (dt == ECFDataType::ENUM_FLAGS)
	{
		int opt = (int)value;
		if (opt < 0) return 0;
		int size = p->metadata.options.size();
		if (opt >= size) return size - 1;
		return opt;
	}
	else if (dt == ECFDataType::BOOL)
	{
		int val = (int)value;
		if (val < 0) return 0;
		if (val > 1) return 1;
		return val;
	}

	return (int)value;
}

double ParameterManager::_checkParameter_no_rounding(int id, double value)
{
	ECFParameter* p = this->m_params[id];
	ECFDataType dt = p->metadata.datatype.front();

	if (dt == ECFDataType::POSITIVE_INTEGER)
	{
		if (value < 1) return 1;
	}
	else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
	{
		if (p->metadata.range)
		{
			if (value < std::atoi(p->metadata.range->min.c_str()))
			{ return std::atoi(p->metadata.range->min.c_str()); }
			else if (value > std::atoi(p->metadata.range->max.c_str()))
			{ return std::atoi(p->metadata.range->max.c_str()); }
		}
		else if (value < 0) return 0;
	}
	else if (dt == ECFDataType::OPTION || dt == ECFDataType::ENUM_FLAGS)
	{
		if (value < 0) return 0;
		int size = p->metadata.options.size();
		if (value >= size) return size - 1;
	}
	else if (dt == ECFDataType::BOOL)
	{
		if (value < 0) return 0;
		if (value > 1) return 1;
	}

	return value;
}

OutputManager::OutputManager(const AutoOptimizationConfiguration& configuration)
: m_config(configuration) {}

void OutputManager::writeConfiguration(const char type,
	std::vector<double>& configuration, int dimension) const
{
	std::stringstream ss;

	switch(info::ecf->output.logger) {
		case OutputConfiguration::LOGGER::USER:
			ss << type << ",";
			for (int d = 0; d < dimension; ++d)
			{ss << configuration[d] << ",";}
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: CONFIGURATION  %36s | -\n",
			ss.str().c_str());
			break;
		case OutputConfiguration::LOGGER::PARSER:
		default:
			ss << "autoopt: configuration=" << type << ",";
			for (int d = 0; d < dimension; ++d)
			{ss << configuration[d] << ",";}
			ss << "\n";
			eslog::info(ss.str().c_str());
	}
}

void OutputManager::writeAlgorithm() const
{
	switch(info::ecf->output.logger) {
		case OutputConfiguration::LOGGER::USER:
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: ALGORITHM                   %23s | -\n",
			m_config.ecfdescription->getParameter(&m_config.algorithm)->getValue().c_str());
			break;
		case OutputConfiguration::LOGGER::PARSER:
		default:
			std::stringstream ss;
			ss << "autoopt: algorithm=";
			ss << m_config.ecfdescription->getParameter(&m_config.algorithm)->getValue();
			ss << "\n";
			eslog::info(ss.str().c_str());
	}
}

void OutputManager::writeEvaluation(double evaluation) const
{
	std::stringstream ss;

	switch(info::ecf->output.logger) {
		case OutputConfiguration::LOGGER::USER:
			ss << evaluation;
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: EVALUATION                   %22s | -\n",
			ss.str().c_str());
			break;
		case OutputConfiguration::LOGGER::PARSER:
		default:
			ss << "autoopt: evaluation=" << evaluation << "\n";
			eslog::info(ss.str().c_str());
	}
}

void OutputManager::writeTime(double time) const
{
	std::stringstream ss;

	switch(info::ecf->output.logger) {
		case OutputConfiguration::LOGGER::USER:
			ss << time;
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: TIME                   %28s | -\n",
			ss.str().c_str());
			break;
		case OutputConfiguration::LOGGER::PARSER:
		default:
			ss << "autoopt: time=" << time << "\n";
			eslog::info(ss.str().c_str());
	}
}

OptimizationProxy::OptimizationProxy(
	std::vector<ECFParameter*>& parameters,
	const AutoOptimizationConfiguration& configuration
)
: m_params(parameters), m_config(configuration),
  dimension(parameters.size()),
  m_manager(m_params, configuration.population, configuration.rounding_immediate),
  m_output(configuration)
{
	if (info::mpi::rank != 0) return;

	switch(m_config.algorithm)
	{
	case AutoOptimizationConfiguration::ALGORITHM::PARTICLE_SWARM:
		this->m_alg = new PSOAlgorithm(m_manager, m_output,
			m_config.population, m_config.particle_swarm.generations,
			m_config.particle_swarm.C1, m_config.particle_swarm.C2,
			m_config.particle_swarm.W_START, m_config.particle_swarm.W_END);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::DIFFERENTIAL_EVOLUTION:
		this->m_alg = new DEAlgorithm(m_manager, m_output,
			m_config.population,
			m_config.differential_evolution.F,
			m_config.differential_evolution.CR);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::SOMAT3A:
		this->m_alg = new SOMAT3AAlgorithm(m_manager, m_output);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::SOMA:
		this->m_alg = new SOMAAlgorithm(m_manager, m_output,
			m_config.population, m_config.soma.PRT,
			m_config.soma.STEP, m_config.soma.PATH_LENGTH);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::RANDOM:
		this->m_alg = new RandomAlgorithm(m_manager, m_output);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::ALL_PERMUTATIONS:
		this->m_alg = new AllPermutationsAlgorithm(m_manager, m_output);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::MICRO_DER:
		this->m_alg = new MicroDERAlgorithm(m_manager, m_output,
			m_config.micro_der.F_START, m_config.micro_der.F_END,
			m_config.micro_der.CR, m_config.micro_der.RESTART_LIMIT);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::IMPROVED_MICRO_PSO:
		this->m_alg = new ImprovedMicroPSOAlgorithm(m_manager, m_output,
			m_config.population, m_config.improved_micro_pso.C1,
			m_config.improved_micro_pso.C2, m_config.improved_micro_pso.W_START,
			m_config.improved_micro_pso.W_END, m_config.improved_micro_pso.POPULATION_CONVERGENCE,
			m_config.improved_micro_pso.CONVERGENCE_THRESHOLD, m_config.improved_micro_pso.M,
			m_config.improved_micro_pso.BETA, m_config.improved_micro_pso.RHO_START,
			m_config.improved_micro_pso.S_C, m_config.improved_micro_pso.F_C);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::MICRO_PSO:
		this->m_alg = new MicroPSOAlgorithm(m_manager, m_output,
			m_config.population, m_config.micro_pso.C1,
			m_config.micro_pso.C2, m_config.micro_pso.W_START,
			m_config.micro_pso.W_END, m_config.micro_pso.POPULATION_CONVERGENCE,
			m_config.micro_pso.CONVERGENCE_THRESHOLD, m_config.micro_pso.M,
			m_config.micro_pso.BETA, m_config.micro_pso.RHO_START,
			m_config.micro_pso.S_C, m_config.micro_pso.F_C);
		break;
	case AutoOptimizationConfiguration::ALGORITHM::MICRO_GA:
		this->m_alg = new MicroGAAlgorithm(m_manager, m_output, m_config.population,
			m_config.micro_ga.K, m_config.micro_ga.CONVERGENCE_THRESHOLD, 
			m_config.micro_ga.POPULATION_CONVERGENCE);
		break;
	default:;
	}
}

OptimizationProxy::~OptimizationProxy()
{
	delete this->m_alg;
}

void OptimizationProxy::setNextConfiguration()
{
	if (info::mpi::rank == 0) {
		m_last_conf = m_alg->getCurrentSpecimen();
		while (this->isForbidden(m_last_conf))
		{
			m_alg->evaluateCurrentSpecimen(std::numeric_limits<double>::max());
			m_last_conf = m_alg->getCurrentSpecimen(); 
		}
	}
	else {
		m_last_conf.clear();
		m_last_conf.resize(this->dimension);
	}

	Communication::broadcast(m_last_conf.data(), this->dimension,
		MPI_DOUBLE, 0);

	int p = 0;
	for (auto it = m_params.begin(); it != m_params.end(); ++it, ++p)
	{
		ECFDataType dt = (*it)->metadata.datatype.front();
		if (dt == ECFDataType::OPTION || dt == ECFDataType::ENUM_FLAGS)
		{
			(*it)->setValue(
				(*it)->metadata.options[(int)m_last_conf[p]].name
			);
		}
		else if (dt == ECFDataType::BOOL)
		{
			(*it)->setValue(
				((int)m_last_conf[p]) == 0 ? "FALSE" : "TRUE"
			);
		}
		else if (dt == ECFDataType::FLOAT)
		{
			(*it)->setValue(std::to_string(m_last_conf[p]));
		}
		else
		{
			(*it)->setValue(std::to_string((int)m_last_conf[p]));
		}
	}
}

void OptimizationProxy::setConfigurationEvaluation(double value, double realtime)
{
	if (info::mpi::rank == 0)
	{ this->m_alg->evaluateCurrentSpecimen(value); }
	this->m_output.writeTime(realtime);
}

void OptimizationProxy::setConfigurationForbidden(double realtime)
{
	if (info::mpi::rank != 0) return;

	this->m_forbiddens.push_back(m_last_conf);
	this->m_alg->evaluateCurrentSpecimen(std::numeric_limits<double>::max());
	this->m_output.writeTime(realtime);
}

bool OptimizationProxy::isForbidden(std::vector<double>& configuration)
{
	bool forbidden = false;
	for (auto f = m_forbiddens.begin(); f != m_forbiddens.end(); ++f)
	{
		size_t i = 0;
		for (; i < m_params.size(); i++)
		{
			if (!m_manager.areParameterValuesEqual(i, configuration[i], (*f)[i]))
			{ break; }
		}
		if (i == m_params.size()) 
		{ 
			forbidden = true; 
			break; 
		}
	}

	return forbidden;
}
