#ifndef SRC_OPTIMIZATION_PROXY_H_
#define SRC_OPTIMIZATION_PROXY_H_

#include <vector>

#include "../config/configuration.h"
#include "evolution.h"

namespace espreso
{

class EvolutionAlgorithm;

class ParameterManager
{

public:
    ParameterManager(std::vector<ECFParameter*>& parameters);

    int count() const;
    std::vector<double> generateConfiguration();
    double checkParameter(int id, double value);

private:
    std::vector<ECFParameter*>& m_params;
};

enum class OptimizationAlgorithm
{
    PARTICLE_SWARM,
    DIFFERENTIAL_EVOLUTION
};

class OptimizationProxy
{

public:
    OptimizationProxy(
        std::vector<ECFParameter*>& parameters, 
        OptimizationAlgorithm algorithm);
    ~OptimizationProxy();

    void setNextConfiguration();
    void setConfigurationEvaluation(double value);

private:
    std::vector<ECFParameter*>& m_params;
    OptimizationAlgorithm m_algorithm;
    EvolutionAlgorithm* m_alg;

    const int dimension;
    ParameterManager m_manager;

    void setAlgorithm();
};

}

#endif /* SRC_OPTIMIZATION_PROXY_H_ */