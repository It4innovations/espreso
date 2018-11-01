#ifndef SRC_OPTIMIZATION_ALGORITHM_H_
#define SRC_OPTIMIZATION_ALGORITHM_H_

#include "problem.h"

namespace espreso
{

class OptimizationAlgorithm
{
public:
    OptimizationAlgorithm(OptimizationProblem& problem)
    : problem(problem) {}
    virtual ~OptimizationAlgorithm() {}
    
    virtual void run() = 0;

protected:
    OptimizationProblem& problem;
};

class PSOAlgorithm : public OptimizationAlgorithm
{
public:
    PSOAlgorithm(OptimizationProblem& problem);

    void run() override;

private:
    const int population;
    const int generations;

    const float C1;
    const float C2;
    float w;
    const float W_START;
    const float W_END;

    std::vector<std::vector<int> > pBest;
    std::vector<std::vector<float> > velocity;
    std::vector<int> gBest;
};

class DEAlgorithm : public OptimizationAlgorithm
{
public:
    DEAlgorithm(OptimizationProblem& problem);

    void run() override;

private:
    const int population;
    const int generations;

    const float F;
    const float CR;

    std::vector<int> best;
};

}

#endif /* SRC_OPTIMIZATION_ALGORITHM_H_ */