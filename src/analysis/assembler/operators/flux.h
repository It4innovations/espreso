
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct OutputFlux: public ActionOperator {
	OutputFlux(int interval, const ParameterData &dND, const ParameterData &temperature, const ParameterData &conductivity, NamedData *gradient)
	: dND(dND, interval),
	  temp(temperature, interval),
	  conductivity(conductivity, interval),
	  flux(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{

	}

	InputParameterIterator dND, temp, conductivity;
	double* flux;

	void operator++()
	{
		++dND; ++temp; ++conductivity;
		flux += info::mesh->dimension;
	}

	void move(int n)
	{
		dND += n; temp += n; conductivity += n;
		flux += n * info::mesh->dimension;
	}
};

template<size_t nodes, size_t gps>
struct OutputFluxIsotropic2D: public OutputFlux {
	using OutputFlux::OutputFlux;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM2NMN1<nodes>(conductivity[gpindex] / gps, dND.data + 2 * nodes * gpindex, temp.data, flux);
		}
	}
};

template<size_t nodes, size_t gps>
struct OutputFluxIsotropic3D: public OutputFlux {
	using OutputFlux::OutputFlux;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM2NMN1<nodes>(conductivity[gpindex] / gps, dND.data + 3 * nodes * gpindex, temp.data, flux);
		}
	}
};

template<size_t nodes, size_t gps>
struct OutputFlux2D: public OutputFlux {
	using OutputFlux::OutputFlux;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM22M2NMN1<nodes>(1. / gps, conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, temp.data, flux);
		}
	}
};

template<size_t nodes, size_t gps>
struct OutputFlux3D: public OutputFlux {
	using OutputFlux::OutputFlux;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDM33M3NMN1<nodes>(1. / gps, conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, temp.data, flux);
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_FLUX_H_ */
