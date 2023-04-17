
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNELS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNELS_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "esinfo/meshinfo.h"
#include "math/simd/simd.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct BasisKernel: ActionOperator {
	const char* name() const { return "Basis"; }

	BasisKernel()
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}
};

struct CopyCoordinatesKernel: ActionOperator {
	const char* name() const { return "CopyCoordinates"; }

	serializededata<esint, esint>::const_iterator enodes, end;
	bool toGPs;

	CopyCoordinatesKernel()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  toGPs(false)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, bool toGPs)
	{
		this->enodes = enodes;
		this->end = end;
		this->toGPs = toGPs;
	}
};

struct TemperatureKernel: ActionOperator {
	const char* name() const { return "Temperature"; }

	serializededata<esint, esint>::const_iterator enodes, end;
	const double * source;
	bool toGPs;

	TemperatureKernel()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  source(nullptr), toGPs(false)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, const double * source, bool toGPs)
	{
		this->enodes = enodes;
		this->end = end;
		this->source = source;
		this->toGPs = toGPs;
		this->isactive = 1;
	}
};

struct IntegrationKernel: ActionOperator {
	const char* name() const { return "Integration"; }

	IntegrationKernel()
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}
};

struct ExternalExpressionKernel: ActionOperator {
	const char* name() const { return evaluator->expression(); }

	ExternalExpressionKernel(Evaluator *evaluator)
	: evaluator(evaluator)
	{
		isconst = evaluator->isConst();
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	Evaluator *evaluator;

	void setTime(double time, int t)
	{
		evaluator->getTime(t) = time;
	}

	void setFrequency(double frequency, int t)
	{
		evaluator->getFrequency(t) = frequency;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNELS_H_ */
