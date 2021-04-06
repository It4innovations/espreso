
#ifndef SRC_AUTOOPT_EVOLUTION_H_
#define SRC_AUTOOPT_EVOLUTION_H_

#include "config/configuration.h"
#include "proxy.h"

#include <vector>

namespace espreso
{

class ParameterManager;
class OutputManager;

class EvolutionAlgorithm
{
public:
	EvolutionAlgorithm(ParameterManager& manager, OutputManager& output)
	: m_manager(manager), m_output(output) {}
	virtual ~EvolutionAlgorithm() {}

	virtual std::vector<double>& getCurrentSpecimen() = 0;
	virtual void evaluateCurrentSpecimen(double value) = 0;

protected:
	std::vector<std::vector<double> > m_specimens;
	ParameterManager& m_manager;
	OutputManager& m_output;
};

class RandomAlgorithm : public EvolutionAlgorithm
{
public:
	RandomAlgorithm(ParameterManager& manager, OutputManager& output);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	std::vector<double> m_last;
};

class AllPermutationsAlgorithm : public EvolutionAlgorithm
{
public:
	AllPermutationsAlgorithm(ParameterManager& manager, OutputManager& output);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	int m_dimension;
	std::vector<double> m_last;
};

class PSOAlgorithm : public EvolutionAlgorithm
{
public:
	PSOAlgorithm(ParameterManager& manager, OutputManager& output,
		int population, int generations, double C1, double C2,
		double W_START, double W_END);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;
	const int generations;

	int generation;
	std::vector<std::vector<double> >::iterator current;
	bool isInitializing;

	const double C1;
	const double C2;
	double w;
	const double W_START;
	const double W_END;

	std::vector<std::vector<double> > pBest;
	std::vector<std::vector<double> > velocity;
	std::vector<double> gBest;
};

class DEAlgorithm : public EvolutionAlgorithm
{
public:
	DEAlgorithm(ParameterManager& manager, OutputManager& output,
		int population, double F, double CR);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;

	int generation;
	bool isInitializing;
	std::vector<std::vector<double> >::iterator current;
	std::vector<std::vector<double> > new_generation;
	std::vector<double> trial_vector;

	const double F;
	const double CR;

	std::vector<double> best;
};

class SOMAAlgorithm : public EvolutionAlgorithm
{
public:
	SOMAAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, double PRT, double STEP, double PATH_LENGTH);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;

	bool isInitializing;
	std::vector<std::vector<double> >::iterator current;
	double PRT;
	double STEP;
	double PATH_LENGTH;

	std::vector<std::vector<double> > journeys;
	std::vector<std::vector<double> >::iterator current_journey;
	int best_journey;

	std::vector<double> best;
};

class SOMAT3AAlgorithm : public EvolutionAlgorithm
{
public:
	SOMAT3AAlgorithm(ParameterManager& manager,
		OutputManager& output, int p_jumps = 10);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;

	int migration;
	bool isInitializing;
	std::vector<std::vector<double> >::iterator current;
	double FEs;
	double PRT;
	double STEP;
	const double STEP_START;
	const double STEP_END;

	std::vector<std::vector<double>* > Ns;
	std::vector<std::vector<double>* >::iterator current_N;

	std::vector<double>* leader;
	std::vector<std::vector<double> > journeys;
	std::vector<std::vector<double> >::iterator current_journey;
	int best_journey;

	const double JUMPS;
	const double FEs_MAX;
	const int M;
	const int N;
	const int K;

	std::vector<double> best;
};

}

#endif /* SRC_AUTOOPT_EVOLUTION_H_ */
