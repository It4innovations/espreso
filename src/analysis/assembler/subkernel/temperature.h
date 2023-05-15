
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_

#include "subkernel.h"

namespace espreso {

struct Temperature: SubKernel {
	const char* name() const { return "Temperature"; }

	serializededata<esint, esint>::const_iterator enodes, end;
	const double * source;
	bool toGPs;

	Temperature()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  source(nullptr), toGPs(false)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION | Assembler::SOLUTION;
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

template <size_t nodes, size_t gps, class Physics>
struct TemperatureKernel: Temperature, Physics {
	TemperatureKernel(const Temperature &base): Temperature(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				element.temp[n][s] = source[enodes->at(n)];
			}
		}
		if (toGPs) {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.gptemp[gp] = zeros();
				for (size_t n = 0; n < nodes; ++n) {
					element.gptemp[gp] = element.gptemp[gp] + load1(element.N[gp][n]) * element.temp[n];
				}
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_ */
