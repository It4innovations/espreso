#include "proxy.h"

#include <climits>

using namespace espreso;

ParameterManager::ParameterManager(std::vector<ECFParameter*>& parameters)
: m_params(parameters)
{ }

int ParameterManager::count() const
{
    return m_params.size();
}

std::vector<double> ParameterManager::generateConfiguration()
{

    std::vector<double> configuration;
    for (auto it = m_params.begin(); it != m_params.end(); ++it)
    {
        ECFDataType dt = (*it)->metadata.datatype.front();
        if (dt == ECFDataType::INTEGER)
        {
            configuration.push_back( (INT_MIN/2) + rand() );
        }
        else if (dt == ECFDataType::POSITIVE_INTEGER)
        {
            configuration.push_back(rand() % RAND_MAX + 1);
        }
        else if (dt == ECFDataType::NONNEGATIVE_INTEGER)
        {
            configuration.push_back(rand());
        }
        else if (dt == ECFDataType::FLOAT)
        {
            configuration.push_back(
                (double)((INT_MIN/2) + rand()) 
                + ((float) rand()) / (float) RAND_MAX
            );
        }
        else if (dt == ECFDataType::OPTION)
        {
            configuration.push_back(
                rand() % (*it)->metadata.options.size()
            );
        }
        else if (dt == ECFDataType::ENUM_FLAGS)
        {
            configuration.push_back(
                rand() % (*it)->metadata.options.size()
            );
        }
        else if (dt == ECFDataType::BOOL)
        {
            configuration.push_back( rand() % 2 );
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
    return this->_checkParameter_rounded(id, value);
    // return this->_checkParameter_bounds(id, value);
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
    OptimizationAlgorithm algorithm
) 
: m_params(parameters), dimension(parameters.size()), m_manager(parameters),
m_alg(NULL), m_algorithm(algorithm)
{
}

OptimizationProxy::~OptimizationProxy()
{
    if (this->m_alg)
    { delete this->m_alg; }
}

void OptimizationProxy::setNextConfiguration()
{
    this->setAlgorithm();
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

void OptimizationProxy::setAlgorithm()
{
    if (this->m_alg) return;

    switch(m_algorithm)
    {
    case OptimizationAlgorithm::PARTICLE_SWARM:
        this->m_alg = new PSOAlgorithm(m_manager);
        break;
    case OptimizationAlgorithm::DIFFERENTIAL_EVOLUTION:
        this->m_alg = new DEAlgorithm(m_manager);
        break;

    default:;
    }
}