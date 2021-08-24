
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

	static const Variable& get();

	static void print();

	Variable(): offset{0}, increment{0}, val{nullptr}, isconst{0}, update{1}, isconst_interval{nullptr}, update_interval{nullptr} {}
	Variable(int offset, int increment, double *val, int isconst, int update): offset{offset}, increment{increment}, val{val}, isconst{isconst}, update{update}, isconst_interval{nullptr}, update_interval{nullptr} {}
	Variable(int offset, int increment, double *val, const std::vector<int> &isconst, const std::vector<int> &update): offset{offset}, increment{increment}, val{val}, isconst{0}, update{1}, isconst_interval{&isconst}, update_interval{&update} {}

	int offset, increment;
	double *val;

	int isconst, update;
	const std::vector<int> *isconst_interval, *update_interval;
};

}

#endif /* SRC_BASIS_EXPRESSION_VARIABLE_H_ */
