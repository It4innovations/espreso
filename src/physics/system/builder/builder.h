
#ifndef SRC_PHYSICS_SYSTEM_BUILDER_BUILDER_H_
#define SRC_PHYSICS_SYSTEM_BUILDER_BUILDER_H_

#include "basis/containers/point.h"

namespace espreso {

class LinearSystem;
class WSMPSystem;
class SuperLUSystem;
class PARDISOSystem;
class MKLPDSSSystem;
class HYPRESystem;
class FETISystem;
class AssemblerData;
class SolverData;

struct Builder {
	enum Request : int {
		NONE        = 0,
		K           = 1 << 0, // Stiffness matrix
		C           = 1 << 1, // Stiffness matrix
		M           = 1 << 2, // Mass matrix (only in assembler)
		R           = 1 << 3, // Residual forces (only in assembler)
		f           = 1 << 4, // Right-hand side
		BC          = 1 << 5, // Boundary DIRICHLET condition

		KCM         = K | C | M,
		RBCf        = R | BC | f
	};

	virtual void init(WSMPSystem &system) =0;
	virtual void init(SuperLUSystem &system) =0;
	virtual void init(PARDISOSystem &system) =0;
	virtual void init(MKLPDSSSystem &system) =0;
	virtual void init(HYPRESystem &system) =0;
	virtual void init(FETISystem &system) =0;

	virtual void buildSystem(WSMPSystem &system) =0;
	virtual void buildSystem(SuperLUSystem &system) =0;
	virtual void buildSystem(PARDISOSystem &system) =0;
	virtual void buildSystem(MKLPDSSSystem &system) =0;
	virtual void buildSystem(HYPRESystem &system) =0;
	virtual void buildSystem(FETISystem &system) =0;

	virtual void updateSolution(WSMPSystem &system) =0;
	virtual void updateSolution(SuperLUSystem &system) =0;
	virtual void updateSolution(PARDISOSystem &system) =0;
	virtual void updateSolution(MKLPDSSSystem &system) =0;
	virtual void updateSolution(HYPRESystem &system) =0;
	virtual void updateSolution(FETISystem &system) =0;

	virtual void reset(Request matrices, WSMPSystem &system);
	virtual void reset(Request matrices, SuperLUSystem &system);
	virtual void reset(Request matrices, PARDISOSystem &system);
	virtual void reset(Request matrices, MKLPDSSSystem &system);
	virtual void reset(Request matrices, HYPRESystem &system);
	virtual void reset(Request matrices, FETISystem &system);

	virtual ~Builder() {}

	Request matrices;

	// TODO: move builder to subclasses
	double timeIntegrationConstantK;
	double timeIntegrationConstantC;
	double timeIntegrationConstantM;

	double stiffnessDamping;
	double massDamping;
	double structuralDampingCoefficient;
	Point rotationAxis;

	int AFTSamples;
	double internalForceReduction;

	bool tangentMatrixCorrection;
	bool prestress;
	bool rayleighDamping;
	bool coriolisDamping;
	bool spinSoftening;

	Builder()
	: matrices(Request::NONE),
	  timeIntegrationConstantK(1),
	  timeIntegrationConstantC(0),
	  timeIntegrationConstantM(0),
	  stiffnessDamping(0),
	  massDamping(0),
	  structuralDampingCoefficient(0),
	  AFTSamples(0),
	  internalForceReduction(1),
	  tangentMatrixCorrection(false),
	  prestress(false),
	  rayleighDamping(false),
	  coriolisDamping(false),
	  spinSoftening(false)
	{}

protected:
	void reset(Request matrices, AssemblerData &assembler, SolverData &solver);
};

inline Builder::Request  operator| (Builder::Request  m1, const Builder::Request &m2)
		{ return static_cast<Builder::Request>(static_cast<int>(m1) | static_cast<int>(m2)); }
inline Builder::Request  operator& (Builder::Request  m1, const Builder::Request &m2)
		{ return static_cast<Builder::Request>(static_cast<int>(m1) & static_cast<int>(m2)); }
inline Builder::Request& operator|=(Builder::Request &m1, const Builder::Request &m2)
		{ m1 =   static_cast<Builder::Request>(static_cast<int>(m1) | static_cast<int>(m2)); return m1; }
inline Builder::Request& operator&=(Builder::Request &m1, const Builder::Request &m2)
		{ m1 =   static_cast<Builder::Request>(static_cast<int>(m1) & static_cast<int>(m2)); return m1; }
}

#endif /* SRC_PHYSICS_SYSTEM_BUILDER_BUILDER_H_ */
