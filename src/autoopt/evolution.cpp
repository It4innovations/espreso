
#include "evolution.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>

using namespace espreso;

RandomAlgorithm::RandomAlgorithm(ParameterManager& manager,
	OutputManager& output)
: EvolutionAlgorithm(manager, output)
{
	m_last = m_manager.generateConfiguration();
}

std::vector<double>& RandomAlgorithm::getCurrentSpecimen()
{
	this->m_output.writeConfiguration('N', m_last, m_manager.count());
	return m_last;
}

void RandomAlgorithm::evaluateCurrentSpecimen(double value)
{
	this->m_output.writeEvaluation(value);
	m_last = m_manager.generateConfiguration();
}

AllPermutationsAlgorithm::AllPermutationsAlgorithm(ParameterManager& manager,
	OutputManager& output)
: EvolutionAlgorithm(manager, output), m_dimension(manager.count())
{
	m_last.resize(m_dimension);
	for (int i = 0; i < m_dimension; i++)
	{ m_last[i] = m_manager.getParameterMin(i); }
}

std::vector<double>& AllPermutationsAlgorithm::getCurrentSpecimen()
{
	m_output.writeConfiguration('N', m_last, m_manager.count());
	return m_last;
}

void AllPermutationsAlgorithm::evaluateCurrentSpecimen(double value)
{
	m_output.writeEvaluation(value);
	
	int carry = 1;
	for (int i = m_dimension - 1; i > -1; i--)
	{
		int p = m_last[i] + carry;
		if (p > m_manager.getParameterMax(i))
		{
			m_last[i] = m_manager.getParameterMin(i);
			carry = 1;
		}
		else {
			m_last[i] = p;
			break;
		}
	}
}


PSOAlgorithm::PSOAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, int generations, double C1, double C2, double W_START,
	double W_END)
: EvolutionAlgorithm(manager, output),
  population(population), dimension(manager.count()),
  generations(generations), generation(0), isInitializing(true),
  C1(C1), C2(C2),
  w(1), W_START(W_START), W_END(W_END)
{
	for (int s = 0; s < population; s++)
	{ m_specimens.push_back(m_manager.generateConfiguration()); }
	this->current = m_specimens.begin();

	this->velocity.resize(this->population);
	for (size_t s = 0; s < m_specimens.size(); s++)
	{
		// VELOCITY GENERATION
		std::vector<double> velocity_vector;
		velocity_vector.resize(dimension);
		for (int d = 0; d < dimension; d++)
		{ velocity_vector[d] = m_manager.generateDecimal(); }
		this->velocity[s] = velocity_vector;
	}
}

std::vector<double>& PSOAlgorithm::getCurrentSpecimen()
{
	if (isInitializing)
	{
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else
	{
		m_output.writeConfiguration('N', *current, dimension);
		return *current;
	}
}

void PSOAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		this->pBest.push_back(*current);
		if (!gBest.size()) gBest = *current;

		// SEARCH OF gBest
		if (this->gBest.back() > value)
		{ this->gBest = *current; }

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			this->w = W_START - (((W_START - W_END) * this->generation) / this->generations);
			this->migrateSpecimen();
		}
	} else {
		int s = std::distance(m_specimens.begin(), current);

		(*current)[this->dimension] = value;
		m_output.writeEvaluation(value);

		if ((*current)[this->dimension] < this->pBest[s].back())
		{
			this->pBest[s] = *current;
		}
		if (this->pBest[s].back() < this->gBest.back())
		{
			this->gBest = this->pBest[s];
		}

		// std::cout << "BEST ";
		// for (int i = 0; i < current->size(); i++) std::cout << gBest[i] << " ";
		// std::cout << std::endl;

		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->generation++;
			this->w = W_START - (((W_START - W_END) * this->generation) / this->generations);
		}

		this->migrateSpecimen();
	}
}

void PSOAlgorithm::migrateSpecimen()
{
	double r1 = m_manager.generateDecimal();
	double r2 = m_manager.generateDecimal();
	int s = std::distance(m_specimens.begin(), current);

	for (int v_i = 0; v_i < this->dimension; v_i++)
	{
		this->velocity[s][v_i] =
			w * velocity[s][v_i]
			+ C1 * r1 * (pBest[s][v_i] - m_specimens[s][v_i])
			+ C2 * r2 * (gBest[v_i] - m_specimens[s][v_i]);

		double new_value = m_manager.checkParameter(
			v_i,
			m_specimens[s][v_i] + velocity[s][v_i]
		);
		m_specimens[s][v_i] = new_value;
	}
}


ImprovedMicroPSOAlgorithm::ImprovedMicroPSOAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, double C1, double C2, double W_START,
	double W_END, double pop_convergence, double convergence_threshold, 
	int M, double BETA, double RHO_START, double S_c, double F_c)
: EvolutionAlgorithm(manager, output),
  population(population), dimension(manager.count()),
  generation(0), isInitializing(true),
  C1(C1), C2(C2),
  w(W_START), W_START(W_START), W_END(W_END),
  POP_CONVERGENCE(pop_convergence), CONVERGENCE_THRESHOLD(convergence_threshold),
  M(M), BETA(BETA), RHO_START(RHO_START), rho(RHO_START), best_successes(0), best_S_c(S_c), 
  best_failures(0), best_F_c(F_c), gBest_index(0)
{
	for (int s = 0; s < population; s++)
	{ m_specimens.push_back(m_manager.generateConfiguration()); }
	this->current = m_specimens.begin();

	this->velocity.resize(this->population);
	for (size_t s = 0; s < m_specimens.size(); s++)
	{
		// VELOCITY GENERATION
		std::vector<double> velocity_vector(dimension, 0.0f);
		for (int d = 0; d < dimension; d++)
		{ velocity_vector[d] = m_manager.generateDecimal(); }
		this->velocity[s] = velocity_vector;
	}

	// INITIALIZE PARAMETER BOUNDARIES
	for (int p = 0; p < dimension; p++)
	{
		param_boundaries.push_back(
			{m_manager.getParameterMin(p), m_manager.getParameterMax(p)}
		);
	}
}

std::vector<double>& ImprovedMicroPSOAlgorithm::getCurrentSpecimen()
{
	if (isInitializing)
	{
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else
	{
		m_output.writeConfiguration('N', *current, dimension);
		return *current;
	}
}

void ImprovedMicroPSOAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		this->pBest.push_back(*current);
		if (!gBest.size()) gBest = *current;

		// SEARCH OF gBest
		if (this->gBest.back() > value)
		{ 
			this->gBest = *current; 
			this->gBest_index = std::distance(m_specimens.begin(), current);
		}

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			if (generation == 0) { total_gBest = gBest; }
			this->migrateSpecimen();
		}
	} else {
		int s = std::distance(m_specimens.begin(), current);

		(*current)[this->dimension] = value;
		m_output.writeEvaluation(value);

		// pBest
		if ((*current)[this->dimension] < this->pBest[s].back())
		{
			this->pBest[s] = *current;
			if (s == gBest_index)
			{ this->best_successes++; this->best_failures = 0; }
		}
		else if (s == gBest_index)
		{ this->best_failures++; this->best_successes = 0; }
		
		// gBest
		if (this->pBest[s].back() < this->gBest.back())
		{
			this->gBest = this->pBest[s];
			if (s != this->gBest_index)
			{ this->best_successes = 0; this->best_failures = 0; }
			this->gBest_index = s;
		}

		bool restart = false;
		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->generation++;

			// RHO UPDATE
			if (best_successes > best_S_c) { rho = 2.0f * rho; }
			else if (best_failures > best_F_c) { rho = 0.5f * rho; }
			
			restart = this->detectRestartAndPerform();
		}

		if (!restart) this->migrateSpecimen();
	}
}

double ImprovedMicroPSOAlgorithm::checkParameterBoundaries(int id, double value)
{
	if (param_boundaries[id][0] > value) { return param_boundaries[id][0]; }
	else if (param_boundaries[id][1] < value) { return param_boundaries[id][1]; }
	else { return value; }
}

void ImprovedMicroPSOAlgorithm::migrateSpecimen()
{
	double r1 = m_manager.generateDecimal();
	double r2 = m_manager.generateDecimal();
	double r3 = m_manager.generateDecimal();
	int s = std::distance(m_specimens.begin(), current);

	// REPULSION COMPUTATION
	std::vector<double> repulsion(this->dimension, 0.0f);
	for (size_t l = 0; l < blacklist.size(); l++)
	{
		std::vector<double> d_li(this->dimension, 0.0f);
		for (int i = 0; i < dimension; i++)
		{ d_li[i] = current->at(i) - blacklist[l][i]; }
		
		double d_li_size = 0.0f;
		for (int i = 0; i < dimension; i++)
		{ d_li_size += d_li[i] * d_li[i]; }
		d_li_size = sqrt(d_li_size);
		double d_li_size_m = pow(d_li_size, M + 1);
		
		for (int i = 0; i < dimension; i++)
		{ d_li[i] = d_li[i] / d_li_size_m; }
		
		for (int i = 0; i < dimension; i++)
		{ repulsion[i] += d_li[i]; }
	}
	for (int i = 0; i < dimension; i++)
	{ repulsion[i] *= dimension; }

	// NEW VELOCITY COMPUTATION AND A PARTICLE MOVE
	for (int v_i = 0; v_i < this->dimension; v_i++)
	{
		// BEST PARTICLE (s==gBest_index) OR A REGULAR PARTICLE
		this->velocity[s][v_i] = (s == gBest_index) ?
			w * velocity[s][v_i]
			+ gBest[v_i] - m_specimens[s][v_i]
			+ rho * (1 - 2 * r3) + repulsion[v_i]
			:
			w * velocity[s][v_i]
			+ C1 * r1 * (pBest[s][v_i] - m_specimens[s][v_i])
			+ C2 * r2 * (gBest[v_i] - m_specimens[s][v_i])
			+ repulsion[v_i];

		double new_value = this->checkParameterBoundaries(
			v_i,
			m_specimens[s][v_i] + velocity[s][v_i]
		);
		m_specimens[s][v_i] = new_value;
	}
}

bool ImprovedMicroPSOAlgorithm::detectRestartAndPerform()
{
	// DETECT RESTART
	int n_converged_specimens = 0;
	for (int s = 0; s < population; s++)
	{
		if (s == gBest_index) continue;

		double sum_diffsquared = 0;
		for (int i = 0; i < dimension; i++)
		{
			double diff = m_specimens[s][i] - gBest[i];
			sum_diffsquared += pow(diff, 2);
		}
		double standarddev = sqrt(sum_diffsquared / (double)dimension);
		if (standarddev < CONVERGENCE_THRESHOLD) n_converged_specimens++;
	}

	double ratio_converged = (double)n_converged_specimens / (double)(population - 1);
	// RESTART?
	if (ratio_converged < POP_CONVERGENCE) return false;

	// IS CONVERGED SOLUTION BETTER?
	bool found_better = (gBest.back() < total_gBest.back()) ? true : false;
	if (!found_better) this->blacklist.push_back(gBest);
	else total_gBest = gBest;

	// UPDATE W
	this->w = found_better ? this->W_START : std::min(this->w * (1.0f + BETA), this->W_END);

	// DYNAMIC SPACE ADJUSTMENT STRATEGY
	double delta = found_better ? this->W_START : std::min(this->w * (1.0f + BETA), this->W_END);
	std::vector<std::vector<double> > new_boundaries;
	for (int p = 0; p < dimension; p++)
	{
		double new_P_min = std::max(gBest[p] - (gBest[p] - param_boundaries[p][0]) * delta,
			m_manager.getParameterMin(p));
		double new_P_max = std::min(gBest[p] + (param_boundaries[p][1] - gBest[p]) * delta,
			m_manager.getParameterMax(p));
		new_boundaries.push_back({new_P_min, new_P_max});
	}

	// SELECTIVE RE-INITIALIZATION STRATEGY
	m_specimens.clear();
	for (int s = 0; s < population; s++)
	{
		// POPULATION
		std::vector<double> specimen(dimension, 0.0f);
		double r4 = m_manager.generateDecimal();
		for (int i = 0; i < dimension; i++)
		{
			double val;
			if (new_boundaries[i][0] < param_boundaries[i][0]
				&& param_boundaries[i][1] < new_boundaries[i][1])
			{
				int coin = m_manager.generateInt(0, 1);
				val = (coin == 0) ?
					new_boundaries[i][0] +
					(param_boundaries[i][0] - new_boundaries[i][0]) * r4
					:
					param_boundaries[i][1] +
					(new_boundaries[i][1] - param_boundaries[i][1]) * r4;
			}
			else if (new_boundaries[i][0] < param_boundaries[i][0])
			{ val = new_boundaries[i][0] +
				(param_boundaries[i][0] - new_boundaries[i][0]) * r4; }
			else if (param_boundaries[i][1] < new_boundaries[i][1])
			{ val = param_boundaries[i][1] +
				(new_boundaries[i][1] - param_boundaries[i][1]) * r4; }
			else
			{ val = m_manager.generateParameter(i); }

			specimen[i] = m_manager.checkParameter(i, val);
		}
		m_specimens.push_back(specimen);

		// VELOCITY
		for (int i = 0; i < dimension; i++)
		{ velocity[s][i] = m_manager.generateDecimal(); }
	}
	this->current = m_specimens.begin();
	this->pBest.clear();
	this->gBest.clear(); this->gBest_index = 0;
	this->isInitializing = true;
	this->rho = RHO_START; best_failures = 0; best_successes = 0;
	this->param_boundaries = new_boundaries;

	return true;
}



DEAlgorithm::DEAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, double F, double CR)
: EvolutionAlgorithm(manager, output), population(population), dimension(manager.count()),
  isInitializing(true), F(F), CR(CR)
{
	for (int s = 0; s < population; s++) {
		m_specimens.push_back(m_manager.generateConfiguration());
	}
	this->current = m_specimens.begin();
}

std::vector<double>& DEAlgorithm::getCurrentSpecimen()
{
	if (isInitializing)
	{
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else
	{		
		this->m_output.writeConfiguration('T', trial_vector, dimension);
		// this->m_output.writeConfiguration('N', new_generation.back());
		return trial_vector;
	}
}

void DEAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		if (!this->best.size()) this->best = *current;
		if (this->best.back() > value) this->best = *current;

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			this->mutateSpecimen();
		}
	}
	else
	{
		trial_vector.push_back(value);
		m_output.writeEvaluation(value);

		// SELECTION
		if (trial_vector.back() < (*current).back())
		{ new_generation.push_back(trial_vector); }
		else
		{ new_generation.push_back(*current); }

		if (this->best.back() > new_generation.back().back())
		{ this->best = new_generation.back(); }

		// std::cout << "Best: " << this->best.back() << std::endl;

		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->m_specimens = new_generation;
			new_generation.clear();
		}

		this->mutateSpecimen();
	}
}

void DEAlgorithm::mutateSpecimen()
{
	// Current individual position in the population
	int s = std::distance(m_specimens.begin(), current);

	// MUTATION
	const int MUTATION_PARENTS = 3;
	std::vector<double> parents[MUTATION_PARENTS];
	int p_idx[MUTATION_PARENTS];
	// Parent #1
	while ( ( p_idx[0] = m_manager.generateSpecimenNumber() ) == s );
	// Parent #2
	do { p_idx[1] = m_manager.generateSpecimenNumber(); } 
	while (s == p_idx[1] || p_idx[0] == p_idx[1]);
	// Parent #3
	do { p_idx[2] = m_manager.generateSpecimenNumber(); } 
	while (s == p_idx[2] || p_idx[0] == p_idx[2] || p_idx[1] == p_idx[2]);

	for (int i = 0; i < MUTATION_PARENTS; i++)
	{ parents[i] = m_specimens[ p_idx[i] ]; }

	std::vector<double> noisy_vector;
	noisy_vector.resize(this->dimension);
	for (int d = 0; d < this->dimension; d++)
	{
		noisy_vector[d] = m_manager.checkParameter(
			d,
			parents[0][d] + F * (parents[1][d] - parents[2][d])
		);
	}

	// CROSSOVER
	int J = m_manager.generateParameterNumber();
	trial_vector.clear(); trial_vector.resize(this->dimension);
	for (int i = 0; i < this->dimension; i++) {
		trial_vector[i] = (m_manager.generateDecimal() < CR) || (i == J)
			? noisy_vector[i] : (*current)[i];
	}
}



MicroDERAlgorithm::MicroDERAlgorithm(ParameterManager& manager, OutputManager& output,
	double F_START, double F_END, double CR, double RESTART_LIMIT)
: EvolutionAlgorithm(manager, output), population(5), dimension(manager.count()),
  generation(0), isInitializing(true), F_START(F_START), F_END(F_END),
  F(manager.generateDouble(F_START, F_END)), CR(CR), LSA_L(2), LSA_k(3), 
  LSA_error(1e-8), LSA_x0(manager.count(), 0.0f), RESTART_LIMIT(RESTART_LIMIT), 
  restarts(0), R(1)
{
	for (int s = 0; s < population; s++) {
		m_specimens.push_back(m_manager.generateConfiguration());
	}
	this->current = m_specimens.begin();
}

std::vector<double>& MicroDERAlgorithm::getCurrentSpecimen()
{
	if (isInitializing)
	{
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else
	{
		if (current->size() == this->dimension)
		{
			this->m_output.writeConfiguration('T', *current, dimension);
			return *current;
		}
		else
		{
			this->m_output.writeConfiguration('R', *ray_vector_iterator, dimension);
			return *ray_vector_iterator;
		}
	}
}

void MicroDERAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			generation++;
			this->sortSpecimens(m_specimens);
			current = m_specimens.begin();
			this->mutateSpecimen();
			this->initLineSearch();
			this->produceRayVectors();
		}
	}
	else
	{
		m_output.writeEvaluation(value);

		// Current individual position in the population
		int s = std::distance(m_specimens.begin(), current);

		// INDIVIDUAL = ELITIST
		if (s < 2)
		{
			ray_vector_iterator->push_back(value);

			if (ray_vector_iterator->back() < ray_vectors[ray_best_position].back())
			{ ray_best_position = std::distance(ray_vectors.begin(), ray_vector_iterator); }

			if (++ray_vector_iterator == ray_vectors.end())
			{
				LSA_delta_r = (ray_best_position == 0) || (ray_best_position == 2 * LSA_k) 
					? 2 * LSA_delta_r : LSA_delta_r / LSA_k;
				LSA_x0 = ray_vectors[ray_best_position];
				this->produceRayVectors();
			}

			if (LSA_delta_r <= LSA_error)
			{
				// SELECTION
				if (LSA_x0.back() < current->back())
				{ new_generation.push_back(LSA_x0); }
				else
				{ new_generation.push_back(*current); }

				current++;
				s = std::distance(m_specimens.begin(), current);
				this->mutateSpecimen();
				this->initLineSearch();
				this->produceRayVectors();
			}
		} // TRACKING + RESTART/RANDOM VECTORS
		else
		{
			if (current->size() == this->dimension)
			{ current->push_back(value); }

			new_generation.push_back(*current);
			current++;
		}

		// SKIP ALREADY EVALUATED INDIVIDUALS
		if (s > 1)
		{
			while ( current != m_specimens.end() && current->size() != this->dimension )
			{ 
				new_generation.push_back(*current);
				current++;
			}
		}

		// END OF A GENERATION
		if (current == m_specimens.end())
		{
			generation++;
			this->m_specimens = this->new_generation;
			this->sortSpecimens(m_specimens);
			this->current = this->m_specimens.begin();

			// RESTART POPULATION CHECK POINT
			if (generation % RESTART_LIMIT == 0)
			{ this->restartPopulation(); }

			this->new_generation.clear();
			F = m_manager.generateDouble(F_START, F_END);
			this->mutateSpecimen();
			this->initLineSearch();
			this->produceRayVectors();
		}
	}
}

void MicroDERAlgorithm::mutateSpecimen()
{
	// Current individual position in the population
	int s = std::distance(m_specimens.begin(), current);

	// MUTATION
	const int MUTATION_PARENTS = 3;
	std::vector<double> parents[MUTATION_PARENTS];
	int p_idx[MUTATION_PARENTS];
	// Parent #1
	while ( ( p_idx[0] = m_manager.generateSpecimenNumber() ) == s );
	// Parent #2
	do { p_idx[1] = m_manager.generateSpecimenNumber(); }
	while (s == p_idx[1] || p_idx[0] == p_idx[1]);
	// Parent #3
	do { p_idx[2] = m_manager.generateSpecimenNumber(); }
	while (s == p_idx[2] || p_idx[0] == p_idx[2] || p_idx[1] == p_idx[2]);

	for (int i = 0; i < MUTATION_PARENTS; i++)
	{ parents[i] = m_specimens[ p_idx[i] ]; }

	std::vector<double> noisy_vector(this->dimension, 0.0f);
	for (int d = 0; d < this->dimension; d++)
	{
		noisy_vector[d] = m_manager.checkParameter(
			d,
			parents[0][d] + F * (parents[1][d] - parents[2][d])
		);
	}

	// CROSSOVER
	int J = m_manager.generateParameterNumber();
	trial_vector.clear(); trial_vector.resize(this->dimension);
	for (int i = 0; i < this->dimension; i++) {
		trial_vector[i] = (m_manager.generateDecimal() <= CR) || (i == J)
			? noisy_vector[i] : (*current)[i];
	}
}

void MicroDERAlgorithm::initLineSearch()
{
	this->LSA_x0.assign(this->dimension, 0.0f);
	this->LSA_delta_r = LSA_L / LSA_k;
}

void MicroDERAlgorithm::produceRayVectors()
{
	this->ray_vectors.clear();
	for (int p = 1; p <= (2 * LSA_k + 1); p++)
	{
		std::vector<double> x_p(this->dimension, 0.0f);
		for (int i = 0; i < this->dimension; i++)
		{
			x_p[i] = m_manager.checkParameter(i,
				LSA_x0[i] + LSA_delta_r * (p - LSA_k - 1) * trial_vector[i]);
		}
		this->ray_vectors.push_back(x_p);
	}
	this->ray_vector_iterator = this->ray_vectors.begin();
	this->ray_best_position = 0;
}

void MicroDERAlgorithm::sortSpecimens(std::vector<std::vector<double> >& specimens)
{
	std::sort(specimens.begin(), specimens.end(), 
		[] (std::vector<double>& s1, std::vector<double>& s2) -> bool 
		{ return s1.back() < s2.back(); }
	);
}

void MicroDERAlgorithm::restartPopulation()
{
	// SORT
	this->sortSpecimens(new_generation);

	// ELITISTS
	this->m_specimens.clear();
	this->m_specimens.push_back(new_generation[0]);
	this->m_specimens.push_back(new_generation[1]);

	// TRACKING VECTORS
	std::vector<double> w1(this->dimension, 0);
	std::vector<double> w2(this->dimension, 0);
	auto fill_w = [&, this] (std::vector<double>& w, std::vector<double>& e) 
	{
		double lambda = m_manager.generateDecimal();
		for (int i = 0; i < dimension; i++)
		{
			w[i] = e[i] * R * lambda;
		}
	};
	fill_w(w1, m_specimens[0]);
	fill_w(w2, m_specimens[1]);
	this->m_specimens.push_back(w1);
	this->m_specimens.push_back(w2);

	// RESTART VECTOR
	this->m_specimens.push_back(this->m_manager.generateConfiguration());

	// SET POINTER TO THE FIRST ELITIST
	this->current = this->m_specimens.begin();

	// ADJUST R
	if ( ( (++restarts) % this->dimension ) == 0 )
	{ this->R = this->R * 0.9f; }
}



DEOldAlgorithm::DEOldAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, double F, double CR)
: EvolutionAlgorithm(manager, output), population(population), dimension(manager.count()),
  isInitializing(true), F(F), CR(CR)
{
	for (int s = 0; s < population; s++) {
		m_specimens.push_back(m_manager.generateVerifiedConfiguration());
	}
	this->current = m_specimens.begin();
}

std::vector<double>& DEOldAlgorithm::getCurrentSpecimen()
{
	if (isInitializing)
	{
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else
	{		
		this->m_output.writeConfiguration('T', trial_vector, dimension);
		// this->m_output.writeConfiguration('N', new_generation.back());
		return trial_vector;
	}
}

void DEOldAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		if (!this->best.size()) this->best = *current;
		if (this->best.back() > value) this->best = *current;

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			this->mutateSpecimen();
		}
	}
	else
	{
		trial_vector.push_back(value);
		m_output.writeEvaluation(value);

		if (trial_vector.back() < (*current).back())
		{ new_generation.push_back(trial_vector); }
		else
		{ new_generation.push_back(*current); }

		if (this->best.back() > new_generation.back().back())
		{ this->best = new_generation.back(); }

		// std::cout << "Best: " << this->best.back() << std::endl;

		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->m_specimens = new_generation;
			new_generation.clear();
		}

		this->mutateSpecimen();
	}
}

void DEOldAlgorithm::mutateSpecimen()
{
	// MUTATION
	const int MUTATION_PARENTS = 3;
	std::vector<double> parents[MUTATION_PARENTS];
	for (int p = 0; p < MUTATION_PARENTS; p++) {
		parents[p] = m_specimens[ m_manager.generateSpecimenNumber() ];
	}
	std::vector<double> noisy_vector;
	noisy_vector.resize(this->dimension);
	for (int d = 0; d < this->dimension; d++)
	{
		noisy_vector[d] = m_manager.checkParameter(
			d,
			parents[0][d] + F * (parents[1][d] - parents[2][d])
		);
	}

	// CROSSOVER
	int random = m_manager.generateParameterNumber();
	trial_vector = noisy_vector;
	trial_vector[random] =
	m_manager.generateDecimal() < CR
	? trial_vector[random] : (*current)[random];
}

SOMAAlgorithm::SOMAAlgorithm(ParameterManager& manager, OutputManager& output,
	int population, double PRT, double STEP, double PATH_LENGTH)
: EvolutionAlgorithm(manager, output), population(population),
  dimension(manager.count()), isInitializing(true),
  PRT(PRT), STEP(STEP), PATH_LENGTH(PATH_LENGTH)
{
	for (int s = 0; s < population; s++)
	{ m_specimens.push_back(m_manager.generateConfiguration()); }
	this->current = m_specimens.begin();
}

std::vector<double>& SOMAAlgorithm::getCurrentSpecimen()
{
	if (isInitializing) {
		m_output.writeConfiguration('I', *current, dimension);
		return *current;
	}
	else {
		m_output.writeConfiguration('J', *current_journey, dimension);
		return *current_journey;
	}
}

void SOMAAlgorithm::evaluateCurrentSpecimen(double value)
{
	auto produceJourneys = [&] ()
	{
		this->journeys.clear();
		const int jumps = PATH_LENGTH / STEP;

		for (int j = 1; j <= jumps; j++)
		{
			double step = j*STEP;
			std::vector<double> PRTVector(dimension, 0);
			for (int d = 0; d < this->dimension; d++)
			{ PRTVector[d] =
				m_manager.generateDecimal() < PRT ?
				1 : 0; }
			std::vector<double> journey(dimension + 1, 0);
			for (int d = 0; d < this->dimension; d++)
			{ journey[d] = current->at(d) +
				(best[d] - current->at(d)) * step * PRTVector[d]; }
			journeys.push_back(journey);
		}
		for (auto j = journeys.begin(); j != journeys.end(); j++)
		{
			for (int d = 0; d < dimension; d++)
			{ (*j)[d] = m_manager.checkParameter(d, (*j)[d]); }
		}
		this->current_journey = journeys.begin();
		this->best_journey = 0;
	};

	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		if (!this->best.size()) this->best = *current;
		if (this->best.back() > value) this->best = *current;

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			produceJourneys();
		}
	}
	else
	{
		(*this->current_journey)[dimension] = value;

		m_output.writeEvaluation(value);

		if (value < journeys[best_journey][dimension])
		{ best_journey = std::distance(journeys.begin(), current_journey); }

		if (++current_journey == journeys.end())
		{
			double best_journey_fitness = journeys[best_journey][dimension];
			if (best_journey_fitness < current->at(dimension))
			{ (*current) = journeys[best_journey]; }
			if (best_journey_fitness < best[dimension])
			{ best = journeys[best_journey]; }

			if (++current == m_specimens.end()) current = m_specimens.begin();

			produceJourneys();
		}
	}
}

SOMAT3AAlgorithm::SOMAT3AAlgorithm(ParameterManager& manager, OutputManager& output,
	int p_jumps)
: EvolutionAlgorithm(manager, output),
  population(10), dimension(manager.count()), migration(0), isInitializing(true),
  FEs(0), PRT(0.05), STEP_START(3.0 / p_jumps), STEP_END(1.5 / p_jumps),
  JUMPS(p_jumps), FEs_MAX(1000), M(5), N(2), K(5)
{
	this->STEP = this->STEP_START;
	for (int s = 0; s < population; s++)
	{ m_specimens.push_back(m_manager.generateConfiguration()); }
	this->current = m_specimens.begin();
}

std::vector<double>& SOMAT3AAlgorithm::getCurrentSpecimen()
{
	if (isInitializing) {
		m_output.writeConfiguration('I', *current, dimension);
		return *current; 
	}
	else {
		m_output.writeConfiguration('J', *current_journey, dimension);
		return *current_journey;
	}
}

void SOMAT3AAlgorithm::evaluateCurrentSpecimen(double value)
{
	auto setupMigration = [&] ()
	{
		this->PRT = 0.05 + 0.90 * (FEs / FEs_MAX);
		this->STEP = STEP_START - STEP_END * (FEs / FEs_MAX);
		std::vector<std::vector<double>* > ms;
		auto compare = [&] (std::vector<double>* a, std::vector<double>* b)
			{ return a->back() < b->back(); };
		for (int m = 0; m < M; m++)
		{ ms.push_back(&m_specimens[m_manager.generateSpecimenNumber()]); }
		std::sort(ms.begin(), ms.end(), compare);
		this->Ns.clear();
		for (int n = 0; n < N; n++)
		{ this->Ns.push_back(ms[n]); }
		current_N = this->Ns.begin();
	};

	auto produceLeader = [&] ()
	{
		leader = &m_specimens[m_manager.generateSpecimenNumber()];
		for (int ki = 1; ki < K; ki++)
		{
			std::vector<double> candidate =
				m_specimens[m_manager.generateSpecimenNumber()];
			leader =
				candidate[dimension] < (*leader)[dimension]
				? &candidate : leader;
		}
		if (leader == *current_N)
			return false;
		else
			return true;
	};

	auto produceJourneys = [&] ()
	{
		this->journeys.clear();
		this->journeys.resize(JUMPS);

		for (int j = 0; j < JUMPS; j++)
		{
			double step = j*STEP;
			std::vector<double> PRTVector(dimension, 0);
			for (int d = 0; d < this->dimension; d++)
			{ PRTVector[d] =
				m_manager.generateDecimal() < PRT ?
				1 : 0; }
			std::vector<double> journey(dimension + 1, 0);
			for (int d = 0; d < this->dimension; d++)
			{ journey[d] = (*current_N)->at(d) +
				((*leader)[d] - (*current_N)->at(d)) * step * PRTVector[d]; }
			journeys[j] = journey;
		}
		for (auto j = journeys.begin(); j != journeys.end(); j++)
		{
			for (int d = 0; d < dimension; d++)
			{ (*j)[d] = m_manager.checkParameter(d, (*j)[d]); }
		}
		this->FEs += JUMPS;
		this->current_journey = journeys.begin();
		this->best_journey = 0;
	};

	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeEvaluation(value);
		if (!this->best.size()) this->best = *current;
		if (this->best.back() > value) this->best = *current;

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
			setupMigration();

			while (!produceLeader())
			{
				if (++current_N == Ns.end())
				{
					this->migration++;
					setupMigration();
				}
			}
			produceJourneys();
		}
	}
	else
	{
		(*this->current_journey)[dimension] = value;

		m_output.writeEvaluation(value);

		if (value < journeys[best_journey][dimension])
		{ best_journey = std::distance(journeys.begin(), current_journey); }
		if (++current_journey == journeys.end())
		{
			double best_journey_fitness = journeys[best_journey][dimension];
			if (best_journey_fitness < (*current_N)->at(dimension))
			{ (*(*current_N)) = journeys[best_journey]; }
			if (best_journey_fitness < best[dimension])
			{ best = journeys[best_journey]; }

			if (++current_N == Ns.end())
			{
				this->migration++;
				setupMigration();
			}

			while (!produceLeader())
			{
				if (++current_N == Ns.end())
				{
					this->migration++;
					setupMigration();
				}
			}
			produceJourneys();
		}
	}
}
