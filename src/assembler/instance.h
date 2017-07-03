
#ifndef SRC_ASSEMBLER_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_H_

#include <cstddef>
#include <vector>
#include <fstream>

namespace espreso {

enum class Property;
class SparseMatrix;
class Solution;
enum class REGULARIZATION;
enum class B0_TYPE;

enum Matrices : int {
	NONE        = 0,
	K           = 1 << 0, // Stiffness matrix
	N           = 1 << 1, // Null spaces (N1, N2)
	M           = 1 << 2, // Mass matrix (only in assembler)
	R           = 1 << 3, // Residual forces (only in assembler)
	f           = 1 << 4, // Right-hand side
	B0          = 1 << 5, // Hybrid Total FETI gluing
	B1          = 1 << 6, // Total FETI gluing
	B1c         = 1 << 7, // Simple B1c
	B1duplicity = 1 << 8, // Lambdas duplicity
	primal      = 1 << 9, // Primal solution
	dual        = 1 << 10  // Dual solution (B1 * Lambdas)
};

struct Instance {

	Instance(size_t domains, const std::vector<int> &neighbours);
	Instance(Instance &other, Matrices &share);
	~Instance();

	void computeKernels(REGULARIZATION regularization);
	void assembleB0(B0_TYPE type);
	//void (*computeKernels)(REGULARIZATION regularization);

	void clear();

	size_t domains;
	std::vector<size_t> &domainDOFCount;
	std::vector<Property> &properties;
	std::vector<int> neighbours;

	std::vector<eslocal> clustersMap;

	std::vector<SparseMatrix> &K, &N1, &N2, &RegMat;
	std::vector<SparseMatrix> &M;
	std::vector<std::vector<double> > &R, &f;

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> &B0;
	std::vector<std::vector<esglobal> > &B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> &B1;
	std::vector<std::vector<esglobal> > &B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<esglobal> > &B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > &B1c, &LB, &B1duplicity;

	std::vector<SparseMatrix> &inequality;
	std::vector<std::vector<double> > &inequalityC;

	// blocks types of B1
	enum CONSTRAINT {
		DIRICHLET,
		EQUALITY_CONSTRAINTS,
		INEQUALITY_CONSTRAINTS,
	};

	std::vector<size_t> &block;


	std::vector<std::vector<double> > primalSolution;
	std::vector<std::vector<double> > dualSolution;

	std::vector<Solution*> solutions;

private:
	void _dummyComputeKernels(REGULARIZATION regularization);

	std::vector<SparseMatrix> _K, _M, _N1, _N2, _RegMat;
	std::vector<std::vector<double> > _R, _f;

	std::vector<SparseMatrix> _B0, _B1, _inequality;
	std::vector<std::vector<double> > _B1c, _LB, _B1duplicity, _inequalityC;
	std::vector<std::vector<esglobal> > _B0subdomainsMap, _B1subdomainsMap, _B1clustersMap;

	std::vector<size_t> _domainDOFCount, _block;
	std::vector<Property> _properties;
};

inline Matrices operator~(Matrices m)
{
	return static_cast<Matrices>(~static_cast<int>(m));
}

inline Matrices operator|(Matrices m1, const Matrices &m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
}

inline Matrices& operator|=(Matrices &m1, const Matrices &m2)
{
	m1 = static_cast<Matrices>(static_cast<int>(m1) | static_cast<int>(m2));
	return m1;
}

inline Matrices operator&(Matrices m1, const Matrices &m2)
{
	return static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
}

inline Matrices& operator&=(Matrices &m1, const Matrices &m2)
{
	m1 = static_cast<Matrices>(static_cast<int>(m1) & static_cast<int>(m2));
	return m1;
}

}


#endif /* SRC_ASSEMBLER_INSTANCE_H_ */
