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
    ParameterManager(std::vector<ECFParameter*>& parameters, bool roundingImmediate = false);

    int count() const;
    int generateInt();
    double generateDouble();
    std::vector<double> generateConfiguration();
    double checkParameter(int id, double value);

private:
    std::vector<ECFParameter*>& m_params;
    const bool m_immediate;
    std::minstd_rand m_generator;
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