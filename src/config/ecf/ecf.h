
#ifndef SRC_CONFIG_ECF_ECF_H_
#define SRC_CONFIG_ECF_ECF_H_

#include "output.h"

#include "pythontestgenerator.h"

#include "input/input.h"
#include "input/generator.h"
#include "input/feti4ilibrary.h"

#include "meshmorphing.h"

#include "physics/physics.h"
#include "physics/coupled.h"
#include "physics/heattransfer.h"
#include "physics/structuralmechanics.h"

#include "config/holders/range.h"

namespace espreso {

struct FunctionDefinition: public ECFDescription {

	enum class AGGREGATOR {
		NONE,
		MIN,
		MAX,
		AVG,
		ABSMIN,
		ABSMAX
	};

	int loadstep;
	std::string region;
	AGGREGATOR aggregator;

	ECFExpression function;

	FunctionDefinition();
};

struct ECF: public ECFDescription {

	enum class INPUT_TYPE {
		EXTERNAL_FILE,
		GENERATOR
	};

	static void init();
	static void init(int *argc, char ***argv, const std::string &app);
	static void init(const std::string &file);
	static void finish();

	PhysicsConfiguration* getPhysics() { return const_cast<PhysicsConfiguration*>(_getPhysics()); }
	const PhysicsConfiguration* getPhysics() const { return _getPhysics(); }

	PythonTestGenerator python_test_generator;

	std::map<size_t, std::string> default_args;
	std::map<std::string, std::string> variables;
	std::map<std::string, ECFRange> ranges;
	std::map<std::string, FunctionDefinition> functions;

	FETI4ILibraryConfiguration feti4ilibrary;

	INPUT_TYPE input_type;
	InputConfiguration input;
	InputGeneratorConfiguration generator;

	MeshMorphing mesh_morphing;

	PhysicsConfiguration::TYPE physics;
	ThermoElasticityConfiguration thermo_elasticity_2d;
	ThermoElasticityConfiguration thermo_elasticity_3d;
	HeatTransferConfiguration heat_transfer_2d;
	HeatTransferConfiguration heat_transfer_3d;
	StructuralMechanicsConfiguration structural_mechanics_2d;
	StructuralMechanicsConfiguration structural_mechanics_3d;

	OutputConfiguration output;

	std::string ecffile;
	std::string exe;
	std::string name;
	std::string outpath;

	ECF();
	ECF(const std::string &file);
	ECF(int *argc, char ***argv, const std::string &app);
	void fill(const std::string &file);
	void fill(int *argc, char ***argv, const std::string &app);
	void set();

protected:
	void _init();

	const PhysicsConfiguration* _getPhysics() const;
};

}

#endif /* SRC_CONFIG_ECF_ECF_H_ */
