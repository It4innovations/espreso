
#ifndef SRC_BASIS_EXPRESSION_VARIABLE_H_
#define SRC_BASIS_EXPRESSION_VARIABLE_H_

#include <string>
#include <vector>
#include <map>

namespace espreso {

struct ECFExpression;

struct Variable {
	static struct {
		struct Region {
			std::map<std::string, Variable> enodes, egps;
		};

		std::map<std::string, Variable> global, element, enodes, egps, node;
		std::vector<Region> region;
	} list;

	static void gather(size_t regions);
	static void analyze(ECFExpression &expr, size_t region);
	static bool create(ECFExpression &expr);
	static bool create(ECFExpression &expr, size_t region);
	static void print();

	Variable(): offset{0}, increment{0}, val{nullptr} {}
	Variable(int offset, int increment, double *val): offset{offset}, increment{increment}, val{val} {}

//	std::string name;

	int offset, increment;
	double *val;
};

}

#endif /* SRC_BASIS_EXPRESSION_VARIABLE_H_ */
