
#ifndef SRC_PHYSICS_LINEARSYSTEM_LINEARSYSTEM_H_
#define SRC_PHYSICS_LINEARSYSTEM_LINEARSYSTEM_H_

#include <cstddef>

namespace espreso {

class Matrix;
class Vectors;
class Composer;
class Builder;

class SystemSolver {
public:
	virtual ~SystemSolver() {}

	virtual void init() =0;
	virtual bool update() =0;
	virtual bool solve() =0;

	virtual double& precision() =0;
};

struct AssemblerData {
	Matrix *K, *M, *C, *CM;
	Vectors *R, *f, *x;
	Vectors *BC;

	AssemblerData(Matrix *K, Matrix *M, Matrix *C, Matrix *CM, Vectors *R, Vectors *f, Vectors *x, Vectors *BC);
	virtual ~AssemblerData();

	virtual void print(const Builder *builder, const char* prefix, const char* suffix) =0;

	Composer *composer;
};

struct SolverData {
	Matrix *K;
	Vectors *R, *f, *x, *y;
	Vectors *BC;

	SolverData(Matrix *K, Vectors *R, Vectors *f, Vectors *x, Vectors *y, Vectors *BC, SystemSolver *linearSolver);
	virtual ~SolverData();

	virtual void setDirichlet(const Builder *builder) =0;

	virtual bool hasReactionForces() { return false; }
	virtual void printData(const Builder *builder, const char* prefix) =0;
	virtual void printSolution(const Builder *builder, const char* prefix) =0;

	SystemSolver *linearSolver;
};

struct LinearSystem {
	Builder *builder;

	virtual int nassemblers() =0;
	virtual int nsolvers() =0;

	virtual AssemblerData* assembler(int index = 0) =0;
	virtual SolverData* solver(int index = 0) =0;

	void init();
	void nextSubstep();
	void assemble();
	void setDirichlet();
	bool solve();
	void solutionChanged();
	void processSolution();

	virtual ~LinearSystem();

protected:
	virtual void _builderInit() =0;
	virtual void _builderReset() =0;
	virtual void _builderCreateSystem() =0;
	virtual void _builderUpdateSolution() =0;

};

}



#endif /* SRC_PHYSICS_LINEARSYSTEM_LINEARSYSTEM_H_ */
