
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct CopyConductivity: public ActionOperator {
	CopyConductivity(int interval, const ParameterData &input, ParameterData &output)
	: input(input, interval),
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

	void move(int n)
	{
		input += n;
		output += n;
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = 0;
			output[ogp + 2] = 0;                   output[ogp + 3] = input[igp + 1];
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = 0;                   output[ogp + 2] = 0;
			output[ogp + 3] = 0;                   output[ogp + 4] = input[igp + 1]; output[ogp + 5] = 0;
			output[ogp + 6] = 0;                   output[ogp + 7] = 0;                   output[ogp + 8] = input[igp + 2];
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = input[igp + 2];
			output[ogp + 2] = input[igp + 2]; output[ogp + 3] = input[igp + 1];
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = input[igp + 3]; output[ogp + 2] = input[igp + 5];
			output[ogp + 3] = input[igp + 3]; output[ogp + 4] = input[igp + 1]; output[ogp + 5] = input[igp + 4];
			output[ogp + 6] = input[igp + 5]; output[ogp + 7] = input[igp + 4]; output[ogp + 8] = input[igp + 2];
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = input[igp + 2];
			output[ogp + 2] = input[igp + 3]; output[ogp + 3] = input[igp + 1];
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
			output[ogp + 0] = input[igp + 0]; output[ogp + 1] = input[igp + 3]; output[ogp + 2] = input[igp + 5];
			output[ogp + 3] = input[igp + 6]; output[ogp + 4] = input[igp + 1]; output[ogp + 5] = input[igp + 4];
			output[ogp + 6] = input[igp + 8]; output[ogp + 7] = input[igp + 7]; output[ogp + 8] = input[igp + 2];
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_ */
