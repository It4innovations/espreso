
#ifndef SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_

#include <vector>
#include <cstddef>

namespace espreso {

enum class Property;
class Instance;
struct Step;

struct InequalityConstraints
{
	static void insertLowerBoundToB1(Instance &instance, const Step &step, const std::vector<Property> &eDOFs, const std::vector<Property> &boundDOFs);

	static void removePositive(Instance &instance, const Step &step, const std::vector<std::vector<double> > &solution, double rho);
	static void reconstruct(Instance &instance, const Step &step);
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_INEQUALITYCONSTRAINTS_H_ */
