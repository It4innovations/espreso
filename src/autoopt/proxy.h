
#ifndef SRC_AUTOOPT_PROXY_H_
#define SRC_AUTOOPT_PROXY_H_

#include "evolution.h"

#include <vector>
#include <random>

namespace espreso
{

class EvolutionAlgorithm;
struct AutoOptimizationConfiguration;

class ParameterManager
{

public:
	ParameterManager(
		std::vector<ECFParameter*>& parameters,
		int population, bool roundingImmediate = false);

	int count() const;
	int generateInt(int start, int endInclusive);
	int generateSpecimenNumber();
	int generateParameterNumber();
	double generateDecimal();
	std::vector<double> generateConfiguration();
	double checkParameter(int id, double value);
	bool areParameterValuesEqual(int id, double val1, double val2);
	double getParameterMin(int id);
	double getParameterMax(int id);

private:
	std::vector<ECFParameter*>& m_params;
	const bool m_immediate;

	std::mt19937 m_generator;
	std::uniform_real_distribution<double> m_dist_decimal;
	std::uniform_int_distribution<int> m_dist_dimension;
	std::uniform_int_distribution<int> m_dist_population;

	double _checkParameter_immediate_rounding(int id, double value);
	double _checkParameter_no_rounding(int id, double value);
};

class OutputManager
{
public:
	OutputManager(const AutoOptimizationConfiguration& configuration);

	void writeConfiguration(const char type, std::vector<double>& configuration);

private:
	const AutoOptimizationConfiguration& m_config;
};

class OptimizationProxy
{

public:
	OptimizationProxy(
		std::vector<ECFParameter*>& parameters,
		const AutoOptimizationConfiguration& configuration);
	~OptimizationProxy();

	void setNextConfiguration();
	void setConfigurationEvaluation(double value);
	void setConfigurationForbidden();

private:
	std::vector<ECFParameter*> m_params;
	const AutoOptimizationConfiguration& m_config;
	EvolutionAlgorithm* m_alg;

	const int dimension;
	ParameterManager m_manager;
	OutputManager m_output;
	std::vector<std::vector<double> > m_forbiddens;

	bool isForbidden(std::vector<double>& configuration);
};

}

#endif /* SRC_AUTOOPT_PROXY_H_ */
