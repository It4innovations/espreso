
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

	void migrateSpecimen();
};

// Source papers:
// [1] microPSO, DOI: 10.1016/j.amc.2006.01.088
// [2] Definition of rho, s_c, f_c, DOI: 10.1109/SIS.2003.1202274
class MicroPSOAlgorithm : public EvolutionAlgorithm
{
public:
	MicroPSOAlgorithm(ParameterManager& manager, OutputManager& output,
		int population, double C1, double C2,
		double W_START, double W_END, double pop_convergence, 
		double convergence_threshold, int M, double BETA, 
		double RHO_START, double S_c, double F_c);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;

	int generation;
	std::vector<std::vector<double> >::iterator current;
	bool isInitializing;

	const double C1;
	const double C2;
	double w;
	const double W_START;
	const double W_END;
	const double POP_CONVERGENCE;
	const double CONVERGENCE_THRESHOLD;
	const int M;
	const double BETA;
	const double RHO_START;
	double rho;
	int best_successes;
	const int best_S_c;
	int best_failures;
	const int best_F_c;

	std::vector<std::vector<double> > pBest;
	std::vector<std::vector<double> > velocity;
	std::vector<double> gBest;
	int gBest_index;
	std::vector<double> total_gBest;
	std::vector<std::vector<double> > blacklist;

	void migrateSpecimen();
	bool detectRestartAndPerform();
};

// Source papers:
// [1] Improved microPSO, DOI: 10.1002/etep.1704
// [2] Original microPSO, DOI: 10.1016/j.amc.2006.01.088
// [3] Definition of rho, s_c, f_c, DOI: 10.1109/SIS.2003.1202274
class ImprovedMicroPSOAlgorithm : public EvolutionAlgorithm
{
public:
	ImprovedMicroPSOAlgorithm(ParameterManager& manager, OutputManager& output,
		int population, double C1, double C2,
		double W_START, double W_END, double pop_convergence, 
		double convergence_threshold, int M, double BETA, 
		double RHO_START, double S_c, double F_c);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;

	int generation;
	std::vector<std::vector<double> >::iterator current;
	bool isInitializing;

	const double C1;
	const double C2;
	double w;
	const double W_START;
	const double W_END;
	const double POP_CONVERGENCE;
	const double CONVERGENCE_THRESHOLD;
	const int M;
	const double BETA;
	const double RHO_START;
	double rho;
	int best_successes;
	const int best_S_c;
	int best_failures;
	const int best_F_c;

	std::vector<std::vector<double> > pBest;
	std::vector<std::vector<double> > velocity;
	std::vector<double> gBest;
	int gBest_index;
	std::vector<double> total_gBest;
	std::vector<std::vector<double> > blacklist;
	std::vector<std::vector<double> > param_boundaries;

	double checkParameterBoundaries(int id, double value);
	void migrateSpecimen();
	bool detectRestartAndPerform();
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

	void mutateSpecimen();
};

// Source papers:
// [1] microDER, DOI: 10.1109/ACCESS.2019.2954296
// [2] Ray-ES, DOI: 10.1007/978-3-642-32937-1_37
class MicroDERAlgorithm : public EvolutionAlgorithm
{
public:
	MicroDERAlgorithm(ParameterManager& manager, OutputManager& output,
		double F_START, double F_END, double CR, double RESTART_LIMIT);

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

	const double F_START;
	const double F_END;
	double F;
	const double CR;

	std::vector<std::vector<double> > ray_vectors;
	std::vector<std::vector<double> >::iterator ray_vector_iterator;
	int ray_best_position;
	const double LSA_L;
	const double LSA_k;
	const double LSA_error;
	std::vector<double> LSA_x0;
	double LSA_delta_r;

	const int RESTART_LIMIT;
	int restarts;
	double R;

	void mutateSpecimen();
	void initLineSearch();
	void produceRayVectors();
	void sortSpecimens(std::vector<std::vector<double> >& specimens);
	void restartPopulation();
};

class DEOldAlgorithm : public EvolutionAlgorithm
{
public:
	DEOldAlgorithm(ParameterManager& manager, OutputManager& output,
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

	void mutateSpecimen();
};

// Sources:
// [1] Simple continuous GA:
//   Simon, Dan. Evolutionary optimization algorithms. ISBN 978-0-470-93741-9.
// [2] Tournament strategy: https://cstheory.stackexchange.com/questions/14758/tournament-selection-in-genetic-algorithms
class GAAlgorithm : public EvolutionAlgorithm
{
public:
	GAAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, int K, double CROSSOVER_RATE, double MUTATION_RATE);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;
	const int K;
	const double CROSSOVER_RATE;
	const double MUTATION_RATE;
	int generation;

	std::vector<std::vector<double> >::iterator current;

	void produceChildren();
};

// Sources:
// [1] microGA, DOI: 10.1117/12.969927 (original paid paper)
// Free alternatives with a short description only
// [2] Golshan, A. et al., A Comparison of Optimization Methods in Cutting
//   Parameters Using Non-dominated Sorting Genetic Algorithm (NSGA-II)
//   and Micro Genetic Algorithm (MGA)
// [3] Coello, C. et al., Multiobjective Optimization using a Micro-Genetic
//   Algorithm, DOI: 10.1007/3-540-44719-9_9
class MicroGAAlgorithm : public EvolutionAlgorithm
{
public:
	MicroGAAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, int K, double CONVERGENCE_THRESHOLD,
	double POP_CONVERGENCE);

	std::vector<double>& getCurrentSpecimen() override;
	void evaluateCurrentSpecimen(double value) override;

private:
	const int population;
	const int dimension;
	const int K;
	const double CONVERGENCE_THRESHOLD;
	const double POP_CONVERGENCE;
	int generation;
	bool isInitializing;

	std::vector<std::vector<double> >::iterator current;

	void produceChildren();
	bool detectRestartAndPerform();
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
