
#ifndef SRC_AUTOOPT_PROXY_H_
#define SRC_AUTOOPT_PROXY_H_

#include "evolution.h"

#include <vector>
#include <random>

namespace espreso
{

class EvolutionAlgorithm;
class AutoOptimizationConfiguration;

class ParameterManager
{

public:
	ParameterManager(
		std::vector<ECFParameter*>& parameters,
		int population, bool roundingImmediate = false);

	int count() const;
	int generateInt(int start, int endInclusive);
	double generateDouble(double start, double endInclusive);
	int generateSpecimenNumber();
	int generateParameterNumber();
	double generateDecimal();
	double generateParameter(int id);
	std::vector<double> generateConfiguration();
	std::vector<double> generateVerifiedConfiguration();
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

	std::vector<double> generateRandomConfiguration();
	std::vector<std::vector<double> > m_verconfs;
};

class OutputManager
{
public:
	OutputManager(const AutoOptimizationConfiguration& configuration);

	void writeConfiguration(const char type, std::vector<double>& configuration, 
		int dimension) const;
	void writeAlgorithm() const;
	void writeEvaluation(double evaluation) const;
	void writeTime(double time) const;

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
	void setConfigurationEvaluation(double value, double realtime);
	void setConfigurationForbidden(double realtime);

private:
	std::vector<ECFParameter*> m_params;
	const AutoOptimizationConfiguration& m_config;
	EvolutionAlgorithm* m_alg;

	const int dimension;
	ParameterManager m_manager;
	OutputManager m_output;
	std::vector<std::vector<double> > m_forbiddens;

	std::vector<double> m_last_conf;
	bool isForbidden(std::vector<double>& configuration);
};

}

#endif /* SRC_AUTOOPT_PROXY_H_ */
