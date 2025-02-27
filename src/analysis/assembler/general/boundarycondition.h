
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_

#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/holders/expression.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"


namespace espreso {

struct BoundaryCondition: public SubKernel {
    // from ecf
    ECFExpression *expression;
    ECFExpressionVector *expressionVector;

    // external
    serializededata<esint, esint>::const_iterator enodes, end;
    double *source;
    esint offset, size;

    BoundaryCondition()
    : expression(nullptr), expressionVector(nullptr),
      enodes(info::mesh->boundaryRegions.front()->nodes->cbegin()), end(info::mesh->boundaryRegions.front()->nodes->cbegin()),
      source(nullptr), offset(0), size(0)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(ECFExpressionVector *expressionVector)
    {
        this->expressionVector = expressionVector;
        if (this->expressionVector) {
            this->isconst = expressionVector->x.evaluator->isConst() && expressionVector->y.evaluator->isConst() && expressionVector->z.evaluator->isConst();
            this->isactive = 1;
            this->needCoordinates = expressionVector->x.evaluator->needCoordinates() && expressionVector->y.evaluator->needCoordinates() && expressionVector->z.evaluator->needCoordinates();
        }
    }

    void activate(ECFExpression *expression)
    {
        this->expression = expression;
        if (this->expression) {
            this->isconst = expression->evaluator->isConst();
            this->isactive = 1;
            this->needCoordinates = expression->evaluator->needCoordinates();
        }
    }

    void activate(double *source, esint offset, esint size, serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end)
    {
        this->source = source;
        this->offset = offset;
        this->size = size;
        this->enodes = enodes;
        this->end = end;
        this->isconst = false;
        this->isactive = 1;
    }

    void activate(double *source, esint offset, esint size)
    {
        this->source = source;
        this->offset = offset;
        this->size = size;
        this->isconst = false;
        this->isactive = 1;
    }
};

struct HarmonicBoundaryCondition: public SubKernel {
    BoundaryCondition magnitude, phase;

    HarmonicBoundaryCondition()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(ECFHarmonicExpressionVector *expressionVector)
    {
        if (expressionVector) {
            magnitude.activate(&expressionVector->magnitude);
            phase.activate(&expressionVector->phase);
            this->isconst = magnitude.isconst && phase.isconst;
            this->isactive = 1;
            this->needCoordinates = magnitude.needCoordinates || phase.needCoordinates;
        }
    }

    void activate(ECFHarmonicExpression *expression)
    {
        if (expression) {
            magnitude.activate(&expression->magnitude);
            phase.activate(&expression->phase);
            this->isconst = magnitude.isconst && phase.isconst;
            this->isactive = 1;
            this->needCoordinates = magnitude.needCoordinates || phase.needCoordinates;
        }
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BOUNDARYCONDITION_H_ */
