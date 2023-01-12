
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct MatrixFiller: public ActionOperator {
	MatrixFiller(int interval, const ParameterData &local, double *global, const esint *position)
	: local(local, interval),
	  global(global), position(position)
	{

	}

	InputParameterIterator local;
	double *global;
	const esint *position;

	void operator++()
	{
		// increment by operator()
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixUpperFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * ((nodes * dimension * (nodes * dimension - 1)) / 2 + nodes * dimension);
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r, local.data += nodes * dimension) {
			for (size_t c = r; c < nodes * dimension; ++c) {
				global[*position++] += *(local.data + c);
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixLowerFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * ((nodes * dimension * (nodes * dimension - 1)) / 2 + nodes * dimension);
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r, local.data += nodes * dimension) {
			for (size_t c = 0; c <= r; ++c) {
				global[*position++] += *(local.data + c);
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct MatrixFullFiller: public MatrixFiller {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * nodes * dimension * nodes * dimension;
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r) {
			for (size_t c = 0; c < nodes * dimension; ++c) {
				global[*position++] += *local.data++;
			}
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct VectorFiller: public ActionOperator {
	VectorFiller(int interval, const ParameterData &rhs, double *global, const esint *position)
	: rhs(rhs, interval),
	  global(global), position(position) {}

	InputParameterIterator rhs;
	double *global;
	const esint *position;

	void operator++()
	{
		++rhs;
	}

	void move(int n)
	{
		rhs += n;
		position += n * nodes * dimension;
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r) {
			global[*position++] += rhs.data[r];
		}
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct VectorSetter: public ActionOperator {
	VectorSetter(int interval, const ParameterData &rhs, double *global, const esint *position, int filter)
	: rhs(rhs, interval),
	  global(global), position(position), filter(filter) {}

	InputParameterIterator rhs;
	int filter;
	double *global;
	const esint *position;

	void operator++()
	{
		++rhs;
	}

	void move(int n)
	{
		rhs += n;
		int dim = 0;
		for (size_t d = 0; d < dimension; ++d) {
			if (filter & (1 << d)) {
				++dim;
			}
		}
		position += n * nodes * dim;
	}

	void operator()()
	{
		for (size_t d = 0; d < dimension; ++d) {
			for (size_t r = 0; r < nodes; ++r) {
				if (filter & (1 << d)) {
					global[*position++] = rhs.data[r];
				}
			}
		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
