
#ifndef SRC_ANALYSIS_ASSEMBLER_CONTROLLER_H_
#define SRC_ANALYSIS_ASSEMBLER_CONTROLLER_H_

#include "operator.h"
#include "parameter.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/expression/variable.h"
#include "mesh/store/nameddata.h"

#include <unordered_map>

namespace espreso {

struct InputHolder {
	virtual int version(int interval) const = 0;
	virtual ~InputHolder() {}
};

struct InputHolderParameterData: public InputHolder {
	const ParameterData &p;

	int version(int interval) const { return p.version[interval]; }
	InputHolderParameterData(const ParameterData &p): p(p) {}
};

//template <typename TEBoundaries, typename TEData>
//struct InputHolderSerializedEData: public InputHolder {
//	const serializededata<TEBoundaries, TEData>* p;
//
//	int version(int interval) const { return 0; } // do we need versioned edata?
//	InputHolderSerializedEData(const serializededata<TEBoundaries, TEData>* p): p(p) {}
//};

struct InputHolderNamedData: public InputHolder {
	const NamedData* p;

	int version(int interval) const { return p->version; }
	InputHolderNamedData(const NamedData* p): p(p) {}
};

struct InputVariable: public InputHolder {
	const Variable* v;

	int version(int interval) const { return v->version_interval ? (*v->version_interval)[interval] : v->update; }
	InputVariable(const Variable* v): v(v) {}
};

struct ParameterStatus {
	ParameterStatus(const ParameterData& parameter): operators(parameter.intervals) {}

	void setConstness(bool constness);

	virtual ~ParameterStatus()
	{
		for (size_t i = 0; i < inputs.size(); ++i) {
			delete inputs[i];
		}
	}

	std::vector<std::vector<ActionOperator*> > operators;
	std::vector<InputHolder*> inputs;
};

class ParameterController {

public:
	void prepare(ParameterData &parameter)
	{
		// calling _get assures that the parameter will be inserted
		auto it = _get(parameter);
		if (it->first->data == nullptr) {
			_get(parameter)->first->resize();
		}
	}

	template <typename ...Other>
	void prepare(ParameterData &parameter, Other&... other)
	{
		prepare(parameter);
		prepare(other...);
	}

	template <typename TEBoundaries, typename TEData>
	void addInput(int interval, ParameterData &parameter, const serializededata<TEBoundaries, TEData>* other)
	{
		_addInput(interval, parameter, other);
	}

	template <typename Input>
	void addInput(ParameterData &parameter, const Input &input)
	{
		_addInput(parameter, input);
	}

	template <typename Input, typename ...Other>
	void addInput(ParameterData &parameter, const Input &input, Other&... other)
	{
		_addInput(parameter, input);
		addInput(parameter, other...);
	}

protected:
	std::unordered_map<ParameterData*, ParameterStatus> parameters;

	std::unordered_map<ParameterData*, ParameterStatus>::iterator _get(ParameterData &parameter)
	{
		auto it = parameters.find(&parameter);
		if (it == parameters.end()) {
			return parameters.insert(std::make_pair(&parameter, ParameterStatus(parameter))).first;
		}
		return it;
	}

	void _addInput(ParameterData &parameter, const ParameterData &other)
	{
		for (size_t i = 0; i < parameter.isconst.size(); ++i) {
			parameter.isconst[i] = parameter.isconst[i] && other.isconst[i];
		}
		_get(parameter)->second.inputs.push_back(new InputHolderParameterData(other));
	}

	void _addInput(ParameterData &parameter, const NamedData* other)
	{
		parameter.setConstness(false);
		_get(parameter)->second.inputs.push_back(new InputHolderNamedData(other));
	}

	template <typename TEBoundaries, typename TEData>
	void _addInput(ParameterData &parameter, const serializededata<TEBoundaries, TEData>* other)
	{
		parameter.setConstness(false);
//		_get(parameter)->second.inputs.push_back(new InputHolderSerializedEData<TEBoundaries, TEData>(other));
	}

	template <typename TEBoundaries, typename TEData>
	void _addInput(int interval, ParameterData &parameter, const serializededata<TEBoundaries, TEData>* other)
	{
		parameter.isconst[interval] = false;
//		_get(parameter)->second.inputs.push_back(new InputHolderSerializedEData<TEBoundaries, TEData>(other));
	}

	void _addInput(ParameterData &parameter, const Variable *v)
	{
//		parameter.isconst[interval] = false;
//		_get(parameter)->second.inputs.push_back(new InputHolderSerializedEData<TEBoundaries, TEData>(other));
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_CONTROLLER_H_ */
