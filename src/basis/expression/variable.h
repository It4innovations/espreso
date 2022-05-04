
#ifndef SRC_BASIS_EXPRESSION_VARIABLE_H_
#define SRC_BASIS_EXPRESSION_VARIABLE_H_

#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "esinfo/stepinfo.h"

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
struct ECFExpression;
struct ECFRange;
struct NamedData;


struct Variable {
	static struct {
		struct Region {
			std::map<std::string, Variable*> enodes, egps;
		};

		std::map<std::string, Variable*> global, element, enodes, egps, node;
		std::vector<Region> region;
	} list;

	static void gather(size_t regions);
	static void analyze(ECFExpression &expr, size_t region);
	static bool create(ECFExpression &expr);
	static bool create(ECFExpression &expr, size_t region);
	static void clear();

	static void print();

	virtual ~Variable() {}

	virtual void set(size_t interval, Evaluator::Params::General &param) const =0;
	virtual int update(size_t interval) const =0;
	virtual int isconst(size_t interval) const =0;

	virtual void updated() =0;
};

struct RangeVariable: public Variable {
	RangeVariable(ECFRange &range);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	double previous;
	ECFRange &range;
};

struct TimeVariable: public Variable {
	TimeVariable(step::Time &time);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	double previous;
	step::Time &time;
};

struct FrequencyVariable: public Variable {
	FrequencyVariable(step::Frequency &frequency);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	double previous;
	step::Frequency &frequency;
};

struct OutputVariable: public Variable {
	OutputVariable(NamedData *data, int offset, int size);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	NamedData *data;
	int offset, size;
};

struct SerializedEdataVariable: public Variable {
	SerializedEdataVariable(serializededata<esint, double> *data, int offset, int size);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	serializededata<esint, double> *data;
	int offset, size;
};

struct SerializedPointsVariable: public Variable {
	SerializedPointsVariable(serializededata<esint, Point> *data, int offset);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	serializededata<esint, Point> *data;
	int offset;
};

struct ParameterVariable: public Variable {
	ParameterVariable(serializededata<esint, double> *data, std::vector<int> &isconst, std::vector<int> &update, int offset, int size);

	void set(size_t interval, Evaluator::Params::General &param) const;
	int update(size_t interval) const;
	int isconst(size_t interval) const;

	void updated();

	serializededata<esint, double> *data;
	std::vector<int> &constness, &updating;
	int offset, size;
};

}

#endif /* SRC_BASIS_EXPRESSION_VARIABLE_H_ */
