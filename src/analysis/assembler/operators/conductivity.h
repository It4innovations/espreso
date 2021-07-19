
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct CopyConductivity: public ActionOperator {
	CopyConductivity(int interval, const ParameterData &input, ParameterData &output)
	: ActionOperator(interval, output.isconst[interval], output.update[interval]),
	  input(input, interval),
	  output(output, interval)
	{

	}

	InputParameterIterator input;
	OutputParameterIterator output;

	void operator++()
	{
		++input;
		++output;
	}

	void reset()
	{

	}
};

template<size_t nodes, size_t gps>
struct CopyDiagonal2DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 2 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = 0;
			this->output.data[ogp + 2] = 0;                         this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	}
};

template<size_t nodes, size_t gps>
struct CopyDiagonal3DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 3 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = 0;                         this->output.data[ogp + 2] = 0;
			this->output.data[ogp + 3] = 0;                         this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = 0;
			this->output.data[ogp + 6] = 0;                         this->output.data[ogp + 7] = 0;                         this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	}
};

template<size_t nodes, size_t gps>
struct CopySymmetric2DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 3 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 2];
			this->output.data[ogp + 2] = this->input.data[igp + 2]; this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	}
};

template<size_t nodes, size_t gps>
struct CopySymmetric3DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 6 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 3]; this->output.data[ogp + 2] = this->input.data[igp + 5];
			this->output.data[ogp + 3] = this->input.data[igp + 3]; this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = this->input.data[igp + 4];
			this->output.data[ogp + 6] = this->input.data[igp + 5]; this->output.data[ogp + 7] = this->input.data[igp + 4]; this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	}
};

template<size_t nodes, size_t gps>
struct CopyAnisotropic2DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 4 * gpindex;
			int ogp = 4 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 2];
			this->output.data[ogp + 2] = this->input.data[igp + 3]; this->output.data[ogp + 3] = this->input.data[igp + 1];
		}
	}
};

template<size_t nodes, size_t gps>
struct CopyAnisotropic3DConductivity: public CopyConductivity {
	using CopyConductivity::CopyConductivity;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			int igp = 9 * gpindex;
			int ogp = 9 * gpindex;
			this->output.data[ogp + 0] = this->input.data[igp + 0]; this->output.data[ogp + 1] = this->input.data[igp + 3]; this->output.data[ogp + 2] = this->input.data[igp + 5];
			this->output.data[ogp + 3] = this->input.data[igp + 6]; this->output.data[ogp + 4] = this->input.data[igp + 1]; this->output.data[ogp + 5] = this->input.data[igp + 4];
			this->output.data[ogp + 6] = this->input.data[igp + 8]; this->output.data[ogp + 7] = this->input.data[igp + 7]; this->output.data[ogp + 8] = this->input.data[igp + 2];
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_ */
