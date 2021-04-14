
#include "evolution.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>

using namespace espreso;

RandomAlgorithm::RandomAlgorithm(ParameterManager& manager,
	OutputManager& output)
: EvolutionAlgorithm(manager, output)
{
	m_last = m_manager.generateConfiguration();
}

std::vector<double>& RandomAlgorithm::getCurrentSpecimen()
{
	this->m_output.writeConfiguration('N', m_last);
	return m_last;
}

void RandomAlgorithm::evaluateCurrentSpecimen(double value)
{
	m_last.push_back(value);
	this->m_output.writeConfiguration('N', m_last);
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
	return m_last;
}

void AllPermutationsAlgorithm::evaluateCurrentSpecimen(double value)
{
	m_last.push_back(value);
	m_output.writeConfiguration('N', m_last);
	m_last.pop_back();
	
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
		return *current;
	}
	else
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

		return *current;
	}
}

void PSOAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeConfiguration('I', *current);
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
		}
	} else {
		int s = std::distance(m_specimens.begin(), current);

		(*current)[this->dimension] = value;

		if ((*current)[this->dimension] < this->pBest[s].back())
		{
			this->pBest[s] = *current;
		}
		if (this->pBest[s].back() < this->gBest.back())
		{
			this->gBest = this->pBest[s];
		}

		m_output.writeConfiguration('N', *current);

		// std::cout << "BEST ";
		// for (int i = 0; i < current->size(); i++) std::cout << gBest[i] << " ";
		// std::cout << std::endl;

		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->generation++;
			this->w = W_START - (((W_START - W_END) * this->generation) / this->generations);
		}
	}
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
		return *current;
	}
	else
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

		return trial_vector;
	}
}

void DEAlgorithm::evaluateCurrentSpecimen(double value)
{
	if (isInitializing)
	{
		current->push_back(value);
		m_output.writeConfiguration('I', *current);
		if (!this->best.size()) this->best = *current;
		if (this->best.back() > value) this->best = *current;

		if (++current == m_specimens.end())
		{
			isInitializing = false;
			current = m_specimens.begin();
		}
	}
	else
	{
		trial_vector.push_back(value);

		if (trial_vector.back() < (*current).back())
		{ new_generation.push_back(trial_vector); }
		else
		{ new_generation.push_back(*current); }

		if (this->best.back() > new_generation.back().back())
		{ this->best = new_generation.back(); }

		this->m_output.writeConfiguration('T', trial_vector);
		this->m_output.writeConfiguration('N', new_generation.back());

		// std::cout << "Best: " << this->best.back() << std::endl;

		if (++current == m_specimens.end())
		{
			current = m_specimens.begin();
			this->m_specimens = new_generation;
			new_generation.clear();
		}
	}
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
	if (isInitializing) return *current;
	else return *current_journey;
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
		m_output.writeConfiguration('I', *current);
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

		m_output.writeConfiguration('J', *current_journey);

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
	if (isInitializing)
	{ return *current; }
	else
	{ return *current_journey; }
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
		m_output.writeConfiguration('I', *current);
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

		m_output.writeConfiguration('J', *current_journey);

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
