
#ifndef SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_

#include <vector>
#include <cstring>

namespace espreso {

class SparseMatrix;
class Mesh;
class ESPRESOSolver;

struct EqualityConstraints;
struct InequalityConstraints;

struct Constraints
{
	friend class EqualityConstraints;
	friend class InequalityConstraints;

	enum BLOCK {
		DIRICHLET = 0,
		EQUALITY_CONSTRAINTS = 1,
		INEQUALITY_CONSTRAINTS = 2,
	};

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;
	std::vector<std::vector<esglobal> > B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<std::vector<esglobal> > B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<esglobal> > B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > B1c;
	std::vector<std::vector<double> > LB;
	std::vector<std::vector<double> > B1duplicity;

	std::vector<SparseMatrix> inequality;
	std::vector<std::vector<double> > inequalityC;

	std::vector<size_t> block;

	std::vector<eslocal> continuityMap;

	Constraints(const ESPRESOSolver &configuration, Mesh &mesh);
	void initMatrices(const std::vector<size_t> &columns);
	void save();



protected:
	static size_t synchronizeOffsets(size_t &offset);

	const ESPRESOSolver &_configuration;
	Mesh &_mesh;
	bool inequalityStored;
};

}


#endif /* SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_ */
