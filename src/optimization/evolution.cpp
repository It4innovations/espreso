#include "evolution.h"

#include <cstdlib>
#include <ctime>

using namespace espreso;

PSOAlgorithm::PSOAlgorithm(ParameterManager& manager) : EvolutionAlgorithm(manager), 
population(10), dimension(manager.count()), generations(9), generation(0), 
C1(2.0f), C2(2.0f), w(1), W_START(0.9f), W_END(0.4f), isInitializing(true)
{
    srand(time(NULL));
    for (int s = 0; s < population; s++)
    { m_specimens.push_back(m_manager.generateConfiguration()); }
    this->current = m_specimens.begin();
    
    this->velocity.resize(this->population);
    for (int s = 0; s < m_specimens.size(); s++)
    {   
        // VELOCITY GENERATION
        std::vector<double> velocity_vector;
        velocity_vector.resize(dimension);
        for (int d = 0; d < dimension; d++)
        { velocity_vector[d] = ((double) rand()) / (double) RAND_MAX; } 
        this->velocity[s] = velocity_vector;
    }
}

std::vector<double> PSOAlgorithm::getCurrentSpecimen()
{
    if (isInitializing)
    {
        return *current;
    }
    else 
    {
        double r1 = ((double) rand()) / (double) RAND_MAX;
        double r2 = ((double) rand()) / (double) RAND_MAX;
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
        std::cout << "I,";
        for (int i = 0; i < current->size(); i++) std::cout << (*current)[i] << ",";
        std::cout << std::endl;
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
    }
    else
    {
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

        std::cout << "M,";
        for (int i = 0; i < current->size(); i++) std::cout << (*current)[i] << ",";
        std::cout << std::endl;

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

DEAlgorithm::DEAlgorithm(ParameterManager& manager) :
EvolutionAlgorithm(manager), population(20), F(0.5f), CR(0.9f),
isInitializing(true), dimension(manager.count())
{
    srand(time(NULL));
    for (int s = 0; s < population; s++)
    { m_specimens.push_back(m_manager.generateConfiguration()); }
    this->current = m_specimens.begin();
}

std::vector<double> DEAlgorithm::getCurrentSpecimen()
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
        for (int p = 0; p < MUTATION_PARENTS; p++)
        { parents[p] = m_specimens[ rand() % this->population ]; }
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
        int random = rand() % this->dimension;
        trial_vector = noisy_vector;
        trial_vector[random] = 
        ( ((double) rand()) / (double) RAND_MAX ) < CR 
        ? trial_vector[random] : (*current)[random];

        return trial_vector;
    }
}

void DEAlgorithm::evaluateCurrentSpecimen(double value)
{
    if (isInitializing)
    {
        current->push_back(value);
        std::cout << "I,";
        for (int i = 0; i < current->size(); i++) std::cout << (*current)[i] << ",";
        std::cout << std::endl;
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

        std::cout << "M,";
        for (int i = 0; i < trial_vector.size(); i++) std::cout << trial_vector[i] << ",";
        std::cout << std::endl;

        // std::cout << "Best: " << this->best.back() << std::endl;

        if (++current == m_specimens.end())
        {
            current = m_specimens.begin();
            this->m_specimens = new_generation;
            new_generation.clear();
        }
    }
}