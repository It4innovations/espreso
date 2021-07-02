
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct Analysis {

template<template<typename> class Solver, typename T, class Analysis, class Configuration>
static AX_LinearSystem<T>* init(Analysis *analysis, Configuration &configuration)
{
	Solver<T> *solver = new Solver<T>(configuration);
	solver->init(analysis);
	analysis->assembler.init();
	analysis->scheme.init(solver);
	analysis->mode.init(solver);
	return solver;
}

};
}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
