#include "proxy.h"

#include <climits>
#include <chrono>

using namespace espreso;

ParameterManager::ParameterManager(std::vector<ECFParameter*>& parameters, 
int population, bool roundingImmediate)
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
//     return static_cast<double>(m_generator()) / static_cast<double>(m_generator.max());
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
            std::cout 
            << "Optimization - UNKNOWN PARAMETER DATA TYPE!" 
            << std::endl;
        }
    }

    return configuration;
}

double ParameterManager::checkParameter(int id, double value)
{
    if (m_immediate) return this->_checkParameter_rounded(id, value);
    else return this->_checkParameter_bounds(id, value);
}

double ParameterManager::_checkParameter_rounded(int id, double value)
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

double ParameterManager::_checkParameter_bounds(int id, double value)
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

OptimizationProxy::OptimizationProxy(
    std::vector<ECFParameter*>& parameters,
    const OptimizationConfiguration& configuration
) 
: m_params(parameters), m_config(configuration), 
dimension(parameters.size()), 
m_manager(parameters, configuration.population, configuration.rounding_immediate)
{
    switch(m_config.algorithm)
    {
    case OptimizationConfiguration::ALGORITHM::PARTICLE_SWARM:
        this->m_alg = new PSOAlgorithm(m_manager, 
            m_config.population, m_config.particle_swarm.generations,
            m_config.particle_swarm.C1, m_config.particle_swarm.C2, 
            m_config.particle_swarm.W_START, m_config.particle_swarm.W_END);
        break;
    case OptimizationConfiguration::ALGORITHM::DIFFERENTIAL_EVOLUTION:
        this->m_alg = new DEAlgorithm(m_manager, m_config.population,
            m_config.differential_evolution.F, 
            m_config.differential_evolution.CR);
        break;
    case OptimizationConfiguration::ALGORITHM::SOMAT3A:
        this->m_alg = new SOMAT3AAlgorithm(m_manager);
        break;
    case OptimizationConfiguration::ALGORITHM::RANDOM:
        this->m_alg = new RandomAlgorithm(m_manager);
    default:;
    }
}

OptimizationProxy::~OptimizationProxy()
{
    delete this->m_alg;
}

void OptimizationProxy::setNextConfiguration()
{
    std::vector<double> configuration = m_alg->getCurrentSpecimen();

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
    this->m_alg->evaluateCurrentSpecimen(value);
}