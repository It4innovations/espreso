#include "problem.h"

#include <cmath>

using namespace espreso;

int OptimizationParameter::minValue() const
{
    return this->m_min;
}

int OptimizationParameter::maxValue() const
{
    return this->m_max;
}

int OptimizationProblem::dimension() const
{
    return this->parameters.size();
}

std::vector<OptimizationParameter>::const_iterator OptimizationProblem::paramBegin() const
{
    return this->parameters.cbegin();
}

std::vector<OptimizationParameter>::const_iterator OptimizationProblem::paramEnd() const
{
    return this->parameters.cend();
}

OptimizationParameter& OptimizationProblem::parameter(int index)
{
    return this->parameters[index];
}

BoothProblem::BoothProblem()
{
    this->parameters.push_back(OptimizationParameter(-10, 10));
    this->parameters.push_back(OptimizationParameter(-10, 10));
}

int BoothProblem::population() const
{
    return 30;
}

int BoothProblem::generations() const
{
    return 10;
}

int BoothProblem::evaluate(const std::vector<int>& specimen)
{
    return pow(specimen[0] + 2 * specimen[1] - 7, 2) + pow(2 * specimen[0] + specimen[1] - 5, 2);
}