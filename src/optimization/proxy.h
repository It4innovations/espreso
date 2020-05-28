#ifndef SRC_OPTIMIZATION_PROXY_H_
#define SRC_OPTIMIZATION_PROXY_H_

#include <vector>
#include <random>

#include "../config/configuration.h"
#include "../config/ecf/solver/optimization/optimization.h"
#include "evolution.h"

namespace espreso
{

class EvolutionAlgorithm;
class OptimizationConfiguration;

class ParameterManager
{

public:
    ParameterManager(std::vector<ECFParameter*>& parameters, 
        int population, bool roundingImmediate = false);

    int count() const;
    int generateInt(int start, int endInclusive);
    int generateSpecimenNumber();
    int generateParameterNumber();
    double generateDecimal();
    std::vector<double> generateConfiguration();
    double checkParameter(int id, double value);

private:
    std::vector<ECFParameter*>& m_params;
    const bool m_immediate;
    
    std::mt19937 m_generator;
    std::uniform_real_distribution<double> m_dist_decimal;
    std::uniform_int_distribution<int> m_dist_dimension;
    std::uniform_int_distribution<int> m_dist_population;

    double _checkParameter_rounded(int id, double value);
    double _checkParameter_bounds(int id, double value);
};

class OptimizationProxy
{

public:
    OptimizationProxy(
        std::vector<ECFParameter*>& parameters, 
        const OptimizationConfiguration& configuration);
    ~OptimizationProxy();

    void setNextConfiguration();
    void setConfigurationEvaluation(double value);

private:
    std::vector<ECFParameter*>& m_params;
    const OptimizationConfiguration& m_config;
    EvolutionAlgorithm* m_alg;

    const int dimension;
    ParameterManager m_manager;

    void setAlgorithm();
};

}

#endif /* SRC_OPTIMIZATION_PROXY_H_ */