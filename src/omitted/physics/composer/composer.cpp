
#include "composer.h"
#include "esinfo/eslog.h"
#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/kernels/kernel.h"
#include "physics/assembler/modules/module.opt.h"
#include "math/matrix.indices.h"
#include "math/matrix.h"
#include "math/vector.sparse.h"
#include "math/vector.dense.h"

#include <cstddef>

using namespace espreso;

Composer::Composer(Kernel *kernel, ModuleOpt *opt)
: kernel(kernel), opt(opt)
{

}

Composer::~Composer()
{
    if (kernel != NULL) { delete kernel; }
}

SolverDataProvider* Composer::provider()
{
    if (kernel != NULL) {
        return kernel->solverDataProvider;
    } else {
        return opt->solverDataProvider;
    }
}

int Composer::solutions()
{
    if (kernel != NULL) {
        return 1;
    } else {
        return 1;
    }
}

void Composer::solutionChanged(Vectors *solution)
{
    double start = eslog::time();
    if (kernel) {
        if (solution) {
            for (esint n = 0; n < solution->nvectors; n++) {
                kernel->solutions[n].fillData(solution->at(n));
            }
        }
        kernel->solutionChanged();
    } else {
        opt->solutionChanged(solution);
    }
    printf("CHANGE SOLUTION: %f\n", eslog::time() - start);
}

void Composer::nextSubstep()
{
    double start = eslog::time();
    if (kernel) {
        kernel->nextSubstep();
    } else {
        opt->nextSubstep();
    }
    printf("ASSEMBLE: %f\n", eslog::time() - start);
}

void Composer::processSolution()
{
    double start = eslog::time();
    if (kernel) {
        kernel->processSolution();
    } else {
        opt->processSolution();
    }
    printf("PROCESS SOLUTION: %f\n", eslog::time() - start);
}

esint Composer::getMatrixSize(esint size, bool omitLower)
{
    return omitLower ? (size * size - size) / 2 + size : size * size;
}

void Composer::insertKPattern(IJ *target, const esint *begin, const esint *end, bool omitLower)
{
    if (omitLower) {
        for (auto row = begin, colbegin = begin; row != end; ++row, ++colbegin) {
            for (auto col = colbegin; col != end; ++col, ++target) {
                if (*row <= *col) {
                    target->row = *row;
                    target->column = *col;
                } else {
                    target->row = *col;
                    target->column = *row;
                }
            }
        }
    } else {
        for (auto row = begin; row != end; ++row) {
            for (auto col = begin; col != end; ++col, ++target) {
                target->row = *row;
                target->column = *col;
            }
        }
    }
}

void Composer::clearMatrices(Builder::Request matrices, AssemblerData *data)
{
    if (matrices & Builder::Request::K) {
        data->K->fill(0);
    }
    if (matrices & Builder::Request::C) {
        data->C->fill(0);
    }
    if (matrices & Builder::Request::M) {
        data->M->fill(0);
    }
    if (matrices & Builder::Request::R) {
        data->R->fill(0);
    }
    if (matrices & Builder::Request::f) {
        data->f->fill(0);
    }
}

void Composer::fillPermutedSparseData(double *target, const std::vector<esint> &indices, const std::vector<esint> &permutation, const std::vector<double> &values)
{
    for (size_t i = 0, j = 0, v = 0; i < indices.size(); i = j, ++v) {
        target[v] = 0;
        while (j < indices.size() && indices[j] == indices[i]) {
            target[v] += values[permutation[j++]];
        }
        target[v] /= j - i;
    }
}

