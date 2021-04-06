
#include "ecf.h"
#include "config/configuration.hpp"
#include "config/reader/reader.h"

#include "basis/utilities/communication.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"

using namespace espreso;

void ECF::init(int *argc, char ***argv)
{
	info::ecf = new ECF(argc, argv);
	if (info::mpi::threading < MPI_THREAD_MULTIPLE) {
		info::ecf->output.mode = OutputConfiguration::MODE::SYNC;
	}
}

void ECF::finish()
{
	delete info::ecf;
}

const PhysicsConfiguration* ECF::_getPhysics() const
{
	switch (physics) {
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D:
		return &thermo_elasticity_2d;
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D:
		return &thermo_elasticity_3d;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		return &heat_transfer_2d;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		return &heat_transfer_3d;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		return &structural_mechanics_2d;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		return &structural_mechanics_3d;
	default:
		eslog::globalerror("Request for unknown physics.");
		return NULL;
	}
}

void ECF::init()
{
	ecfdescription->name = "root";

	REGISTER(python_test_generator, ECFMetaData()
			.setdescription({ "Description of Python test generator (run python tests/generate.py PATH)." }));

	REGISTER(default_args, ECFMetaData()
		.setdescription({ "The index of the argument.", "The argument value." })
		.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::STRING })
		.setpattern({ "0", "VALUE" }));

	REGISTER(variables, ECFMetaData()
			.setdescription({ "A name of variable usable in *.ecf file.", "A value of the variable." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::STRING })
			.setpattern({ "MY_VARIABLE", "VALUE" }));

	ecfdescription->addSpace();

	REGISTER(feti4ilibrary, ECFMetaData()
			.setdescription({ "Settings for FETI4I library." }));

	ecfdescription->addSpace();

	input_type = INPUT_TYPE::EXTERNAL_FILE;
	REGISTER(input_type, ECFMetaData()
			.setdescription({ "Input type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("EXTERNAL_FILE").setdescription("External database."))
			.addoption(ECFOption().setname("GENERATOR").setdescription("Internal mesh generator.")));

	REGISTER(input, ECFMetaData()
			.setdescription({ "Description of Ansys WorkBench format." })
			.allowonly([&] () { return input_type == INPUT_TYPE::EXTERNAL_FILE; }));
	REGISTER(generator, ECFMetaData()
			.setdescription({ "Description of ESPRESO generator." })
			.allowonly([&] () { return input_type == INPUT_TYPE::GENERATOR; }));

	REGISTER(mesh_morphing, ECFMetaData()
			.setdescription({ "Settings for mesh morphing." }));

	physics = PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D;
	REGISTER(physics, ECFMetaData()
			.setdescription({ "Physics" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("THERMO_ELASTICITY_2D").setdescription("Coupled 2D physics."))
			.addoption(ECFOption().setname("THERMO_ELASTICITY_3D").setdescription("Coupled 3D physics."))
			.addoption(ECFOption().setname("HEAT_TRANSFER_2D").setdescription("Heat transfer 2D."))
			.addoption(ECFOption().setname("HEAT_TRANSFER_3D").setdescription("Heat transfer 3D."))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS_2D").setdescription("Structural mechanics 2D."))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS_3D").setdescription("Structural mechanics 3D.")));

	REGISTER(thermo_elasticity_2d, ECFMetaData()
			.setdescription({ "Coupled physics" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D; }));
	REGISTER(thermo_elasticity_3d, ECFMetaData()
			.setdescription({ "Coupled physics" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D; }));
	REGISTER(heat_transfer_2d, ECFMetaData()
			.setdescription({ "Heat transfer 2D" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D; }));
	REGISTER(heat_transfer_3d, ECFMetaData()
			.setdescription({ "Heat transfer 3D" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D; }));
	REGISTER(structural_mechanics_2d, ECFMetaData()
			.setdescription({ "Structural mechanics 2D" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D; }));
	REGISTER(structural_mechanics_3d, ECFMetaData()
			.setdescription({ "Structural mechanics 3D" })
			.allowonly([&] () { return physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }));

	REGISTER(output, ECFMetaData()
			.setdescription({ "Output configurations." }));
}

ECF::ECF()
: mesh_morphing(this),
  thermo_elasticity_2d(DIMENSION::D2),
  thermo_elasticity_3d(DIMENSION::D3),
  heat_transfer_2d(DIMENSION::D2),
  heat_transfer_3d(DIMENSION::D3),
  structural_mechanics_2d(DIMENSION::D2),
  structural_mechanics_3d(DIMENSION::D3),
  output(this)
{
	if (info::ecf == NULL) {
		info::ecf = this;
	} else {
		eslog::globalerror("ESPRESO internal error: cannot create more ECF instances.");
	}
	init();
}

ECF::ECF(const std::string &file)
: ECF()
{
	fill(file);
}

ECF::ECF(int *argc, char ***argv)
: ECF()
{
	fill(argc, argv);
}

void ECF::fill(const std::string &file)
{
	ecffile = file;
	ECFReader::read(*this->ecfdescription, file, this->default_args, this->variables);
}

void ECF::fill(int *argc, char ***argv)
{
	ecffile = ECFReader::read(*this->ecfdescription, argc, argv, this->default_args, this->variables);
}



