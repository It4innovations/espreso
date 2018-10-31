#ifndef SRC_OPTIMIZATION_PROBLEM_H_
#define SRC_OPTIMIZATION_PROBLEM_H_

#include <vector>

namespace espreso
{
    class OptimizationParameter
    {
    public:
        OptimizationParameter(int min_value, int max_value) 
        : m_min(min_value), m_max(max_value) {}
        int minValue() const;
        int maxValue() const;

    private:
        const int m_min;
        const int m_max;
    };

    class OptimizationProblem
    {
    public:
        OptimizationProblem() {}
        virtual ~OptimizationProblem() {}

        virtual int population() const = 0;
        virtual int generations() const = 0;
        int dimension() const;
        std::vector<OptimizationParameter>::const_iterator 
        paramBegin() const;
        std::vector<OptimizationParameter>::const_iterator 
        paramEnd() const;
        OptimizationParameter& parameter(int index);
        virtual int evaluate(const std::vector<int>& specimen) = 0;
    
    protected:
        std::vector<OptimizationParameter> parameters;
    };

    class BoothProblem : public OptimizationProblem
    {
    public:
        BoothProblem();

        int population() const override;
        int generations() const override;
        int evaluate(const std::vector<int>& specimen) override;
    };
}

#endif /* SRC_OPTIMIZATION_PROBLEM_H_ */