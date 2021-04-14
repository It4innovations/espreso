
#include "proxy.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/linearsolver/autoopt.h"
#include "config/ecf/output.h"
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
{ }

int ParameterManager::count() const
{
	return m_params.size();
}

int ParameterManager::generateInt(int start, int endInclusive)
{
	std::uniform_int_distribution<int> d(start, endInclusive);
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
			configuration.push_back(m_generator());
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
		return 0;
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
		return std::numeric_limits<int>::max();
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
		if (val < 0) return 0;
		return val;
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
	else
	{
		return value;
	}
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
		if (value < 0) return 0;
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

const void OutputManager::writeConfiguration(const char type,
	std::vector<double>& configuration)
{
	std::stringstream ss;

	switch(info::ecf->output.logger) {
		case OutputConfiguration::LOGGER::USER:
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: ALGORITHM                   %23s | -\n",
			m_config.ecfdescription->getParameter(&m_config.algorithm)->getValue().c_str());
			ss << type << ",";
			for (auto p = configuration.cbegin();
				p != configuration.cend();
			 	++p)
			{ss << *p << ",";}
			eslog::solver("     - | AUTOMATIC OPTIMIZATION :: CONFIGURATION  %36s | -\n",
			ss.str().c_str());
			break;
		case OutputConfiguration::LOGGER::PARSER:
		default:
			std::stringstream ss;
			ss << "autoopt: algorithm=";
			ss << m_config.ecfdescription->getParameter(&m_config.algorithm)->getValue();
			ss << " configuration=" << type << ",";
			for (auto p = configuration.cbegin();
				p != configuration.cend();
				++p)
			{ss << *p << ",";}
			ss << "\n";
			eslog::info(" = ====================== AUTOMATIC OPTIMIZATION OF SOLVER PARAMETERS ====================== = \n");
			eslog::info(ss.str().c_str());
			eslog::info(" = ========================================================================================= = \n");
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
	default:;
	}
}

OptimizationProxy::~OptimizationProxy()
{
	delete this->m_alg;
}

void OptimizationProxy::setNextConfiguration()
{
	std::vector<double> configuration;

	if (info::mpi::rank == 0) {
		configuration = m_alg->getCurrentSpecimen();
		while (this->isForbidden(configuration))
		{
			m_alg->evaluateCurrentSpecimen(std::numeric_limits<double>::max());
			configuration = m_alg->getCurrentSpecimen(); 
		}
	}
	else {
		configuration.resize(this->dimension);
	}

	Communication::broadcast(configuration.data(), this->dimension,
		MPI_DOUBLE, 0);

	int p = 0;
	for (auto it = m_params.begin(); it != m_params.end(); ++it, ++p)
	{
		ECFDataType dt = (*it)->metadata.datatype.front();
		if (dt == ECFDataType::OPTION || dt == ECFDataType::ENUM_FLAGS)
		{
			(*it)->setValue(
				(*it)->metadata.options[(int)configuration[p]].name
			);
		}
		else if (dt == ECFDataType::BOOL)
		{
			(*it)->setValue(
				((int)configuration[p]) == 0 ? "FALSE" : "TRUE"
			);
		}
		else if (dt == ECFDataType::FLOAT)
		{
			(*it)->setValue(std::to_string(configuration[p]));
		}
		else
		{
			(*it)->setValue(std::to_string((int)configuration[p]));
		}
	}
}

void OptimizationProxy::setConfigurationEvaluation(double value)
{
	double maxtime = 0;
	Communication::reduce(&value, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0);

	if (info::mpi::rank == 0)
	{ this->m_alg->evaluateCurrentSpecimen(maxtime); }
}

void OptimizationProxy::setConfigurationForbidden()
{
	if (info::mpi::rank != 0) return;

	this->m_forbiddens.push_back(m_alg->getCurrentSpecimen());
	this->m_alg->evaluateCurrentSpecimen(std::numeric_limits<double>::max());
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
