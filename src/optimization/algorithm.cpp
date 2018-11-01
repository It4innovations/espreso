#include "algorithm.h"
#include "population.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
 
using namespace espreso;

PSOAlgorithm::PSOAlgorithm(OptimizationProblem& problem)
: OptimizationAlgorithm(problem), 
population(problem.population()), generations(problem.generations()),
C1(0.1f), C2(0.2f), w(1), W_START(0.9f), W_END(0.4f)
{
    srand(time(NULL));
}

void PSOAlgorithm::run()
{
    // INITIALIZATION
    auto specimens = PopulationGenerator::generate(this->problem);
    for (auto s = specimens.begin(); s != specimens.end(); ++s)
    {
        s->push_back(this->problem.evaluate(*s));
    }
    this->pBest = specimens;
    this->gBest = specimens[0];
    this->velocity.resize(this->population);
    for (int s = 0; s < specimens.size(); s++)
    {
        // SEARCH OF gBest
        if (this->gBest.back() > specimens[s].back())
        { this->gBest = specimens[s]; }
        
        // VELOCITY GENERATION
        std::vector<float> velocity_vector;
        velocity_vector.resize(this->problem.dimension());
        for (int d = 0; d < this->problem.dimension(); d++)
        { velocity_vector[d] = ((float) rand()) / (float) RAND_MAX; } 
        this->velocity[s] = velocity_vector;
    }

    // ALGORITHM
    for (int g = 0; g < this->generations; g++)
    {
        this->w = W_START - (((W_START - W_END) * g) / this->generations);
        
        std::cout << "GENERATION " << g + 1 << std::endl;

        for (int s = 0; s < specimens.size(); ++s)
        {
            float r1 = ((float) rand()) / (float) RAND_MAX;
            float r2 = ((float) rand()) / (float) RAND_MAX;

            std::cout << "UNIT " << s << " ";

            for (int v_i = 0; v_i < this->problem.dimension(); v_i++)
            {
                this->velocity[s][v_i] = 
                    w * velocity[s][v_i] 
                    + C1 * r1 * (pBest[s][v_i] - specimens[s][v_i])
                    + C2 * r2 * (gBest[v_i] - specimens[s][v_i]);
                
                int new_value = specimens[s][v_i] + velocity[s][v_i];
                if (new_value < this->problem.parameter(v_i).minValue()) 
                    new_value = this->problem.parameter(v_i).minValue();
                if (new_value > this->problem.parameter(v_i).maxValue())
                    new_value = this->problem.parameter(v_i).maxValue();
                specimens[s][v_i] = new_value;
                
                std::cout << new_value << ", ";
            }

            specimens[s][this->problem.dimension()] = this->problem.evaluate(specimens[s]);

            std::cout << "EVAL: " << specimens[s].back() << std::endl;
           
            if (specimens[s][this->problem.dimension()] < this->pBest[s].back())
            {
                this->pBest[s] = specimens[s];
            }
            if (this->pBest[s].back() < this->gBest.back())
            {
                this->gBest = this->pBest[s];
            }
        }
    }

}

DEAlgorithm::DEAlgorithm(OptimizationProblem& problem)
: OptimizationAlgorithm(problem), 
population(problem.population()), generations(problem.generations()),
F(0.5f), CR(0.9f)
{
    srand(time(NULL));
}

void DEAlgorithm::run()
{
    // INITIALIZATION
    auto specimens = PopulationGenerator::generate(this->problem);
    for (auto s = specimens.begin(); s != specimens.end(); ++s)
    {
        s->push_back(this->problem.evaluate(*s));
    }

    // ALGORITHM
    for (int g = 0; g < this->generations; g++)
    {
        std::cout << "GENERATION " << g << std::endl;

        std::vector<std::vector<int> > new_generation;
        new_generation.resize(this->problem.population());

        for (int s = 0; s < specimens.size(); s++)
        {
            // MUTATION
            const int MUTATION_PARENTS = 3;
            std::vector<int> parents[MUTATION_PARENTS];
            for (int p = 0; p < MUTATION_PARENTS; p++)
            { parents[p] = specimens[ rand() % this->population ]; }
            std::vector<int> noisy_vector;
            noisy_vector.resize(this->problem.dimension());
            for (int d = 0; d < this->problem.dimension(); d++)
            {
                noisy_vector[d] = parents[0][d] + F * (parents[1][d] - parents[2][d]);
                if (noisy_vector[d] < this->problem.parameter(d).minValue())
                { noisy_vector[d] = this->problem.parameter(d).minValue(); }
                if (noisy_vector[d] > this->problem.parameter(d).maxValue())
                { noisy_vector[d] = this->problem.parameter(d).maxValue(); }
            }

            // CROSSOVER
            int random = rand() % this->problem.dimension();
            std::vector<int> trial_vector = noisy_vector;
            trial_vector[random] = 
            ( ((float) rand()) / (float) RAND_MAX ) < CR 
            ? trial_vector[random] : specimens[s][random];
            trial_vector.push_back(this->problem.evaluate(trial_vector));

            if (trial_vector.back() < specimens[s].back())
            { new_generation[s] = trial_vector; }
            else
            { new_generation[s] = specimens[s]; }

            std::cout << "UNIT ";
            for (int d = 0; d < this->problem.dimension(); d++)
            { std::cout << new_generation[s][d] << ", "; }
            std::cout << " -- " << new_generation[s].back() << std::endl;
        }

        specimens = std::move(new_generation);
    }
}