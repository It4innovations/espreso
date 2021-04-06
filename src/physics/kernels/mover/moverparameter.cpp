
#include "moverparameter.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template<>
void MoverParameter::setDataPointer<ElementNodeValues>(ElementNodeValues &parameter, esint offset)
{
	if (parameter.values) {
		if (parameter.isConst) {
			data = parameter.values->datatarray().begin(omp_get_thread_num());
			increment = MoverIncrement();
		} else {
			esint _offset = (parameter.nodes->begin() + offset)->begin() - parameter.nodes->begin()->begin();
			data = (parameter.values->begin() + _offset)->begin();
		}
	}
}

template <>
void MoverParameter::setDataPointer<OutputNodes>(OutputNodes &parameter, esint offset)
{
	if (parameter.data) {
		data = parameter.data->data.data() + offset * parameter.data->dimension;
	}
}

template <>
void MoverParameter::setDataPointer<OutputElements>(OutputElements &parameter, esint offset)
{
	if (parameter.data) {
		data = parameter.data->data.data() + offset * parameter.data->dimension;
	}
}

template<>
void MoverParameter::setOutput<OutputNodes>(
		OutputNodes &output,
		NamedData::DataType datatype,
		const char* name,
		bool guard,
		step::TYPE restriction)
{
	if (guard) {
		output.data = info::mesh->nodes->appendData(dimension, datatype, name, restriction);
		std::fill(output.data->data.begin(), output.data->data.end(), defaultValue);
	}
}

template<>
void MoverParameter::setOutput<OutputElements>(
		OutputElements &output,
		NamedData::DataType datatype,
		const char* name,
		bool guard,
		step::TYPE restriction)
{
	if (guard) {
		output.data = info::mesh->elements->appendData(dimension, datatype, name, restriction);
		std::fill(output.data->data.begin(), output.data->data.end(), defaultValue);
	}
}

void MoverCoordinates::clear()
{
	if (kernel.values) { delete kernel.values; }
}

void MoverCoordinates::initKernelIterator(esint offset)
{
	MoverParameter::setDataPointer(kernel, offset);
}

void MoverCoordinates::initOutputIterator(esint offset)
{
	data = NULL;
}

void MoverCoordinates::updatingFromInput(MoverInstance *iterator)
{
	iterator->now(input, kernel);
}

void MoverCoordinates::updatingFromOutput(MoverInstance *iterator)
{

}

void MoverCoordinates::set(serializededata<esint, Point>* points, serializededata<esint, esint> *nodes)
{
	input.dimension = dimension;
	input.values = points;

	kernel.dimension = dimension;
	kernel.nodes = nodes;
	kernel.values = new serializededata<esint, double>(dimension, nodes->datatarray().distribution());
}

template<>
void MoverOutputParameter<OutputElements>::initKernelIterator(esint offset)
{
	MoverParameter::setDataPointer(output, offset);
}

template<>
void MoverOutputParameter<OutputNodes>::initKernelIterator(esint offset)
{
	data = NULL;
}

template<>
void MoverOutputParameter<OutputElements>::initOutputIterator(esint offset)
{
	data = NULL;
}

template<>
void MoverOutputParameter<OutputNodes>::initOutputIterator(esint offset)
{
	MoverParameter::setDataPointer(output, offset);
}

}
