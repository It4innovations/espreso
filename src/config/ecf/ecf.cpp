
#include "ecf.h"
#include "config/configuration.hpp"
#include "config/reader/reader.h"

#include "basis/logging/timelogger.h"
#include "basis/utilities/sysutils.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/systeminfo.h"

using namespace espreso;

void ECF::init()
{
    info::ecf = new ECF();
}

void ECF::init(int *argc, char ***argv, const std::string &app)
{
    info::ecf = new ECF(argc, argv, app);
    if (info::mpi::threading < MPI_THREAD_MULTIPLE) {
        info::ecf->output.mode = OutputConfiguration::MODE::SYNC;
    }
}

void ECF::init(const std::string &file)
{
    info::ecf = new ECF(file);
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
//    case PhysicsConfiguration::TYPE::THERMO_ELASTICITY:
//        return &thermo_elasticity_2d;
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        return &heat_transfer;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        return &structural_mechanics;
//    case PhysicsConfiguration::TYPE::ACOUSTICS:
//        return &acoustics;
    default:
        eslog::globalerror("Request for unknown physics.");
        return NULL;
    }
}

//FunctionDefinition::FunctionDefinition()
//{
//    loadstep = -1;
//    REGISTER(loadstep, ECFMetaData()
//            .setdescription({ "Expression is evaluated only for a given loadstep (use -1 for the last)." })
//            .setdatatype({ ECFDataType::INTEGER }));
//
//    REGISTER(region, ECFMetaData()
//            .setdescription({ "Expression is evaluated only for a given loadstep (use -1 for the last)." })
//            .setdatatype({ ECFDataType::REGION }));
//
//    aggregator = AGGREGATOR::NONE;
//    REGISTER(aggregator, ECFMetaData()
//            .setdescription({ "Variant" })
//            .setdatatype({ ECFDataType::OPTION })
//            .addoption(ECFOption().setname("NONE").setdescription("Define a new value for all nodes/elements."))
//            .addoption(ECFOption().setname("MIN").setdescription("Minimum across a given region."))
//            .addoption(ECFOption().setname("MAX").setdescription("Maximum across a given region."))
//            .addoption(ECFOption().setname("AVG").setdescription("Average across a given region."))
//            .addoption(ECFOption().setname("ABSMIN").setdescription("Absolute value of the minimun across a given region."))
//            .addoption(ECFOption().setname("ABSMAX").setdescription("Absolute value of the maximun across a given region.")));
//
//    REGISTER(function, ECFMetaData()
//            .setdescription({ "The function definition." })
//            .setdatatype({ ECFDataType::EXPRESSION}));
//}

void ECF::_init()
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

    REGISTER(ranges, ECFMetaData()
            .setdescription({ "A name of range usable in *.ecf file.", "A range of the variable" })
            .setdatatype({ ECFDataType::STRING, ECFDataType::RANGE })
            .setpattern({ "MY_VARIABLE", "MIN : MAX: STEP" }));

//    REGISTER(functions, ECFMetaData()
//            .setdescription({ "A name of a function usable in expressions.", "The function definition." })
//            .setdatatype({ ECFDataType::STRING })
//            .setpattern({ "FNC" }));

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
    REGISTER(coupling, ECFMetaData()
            .setdescription({ "Settings for coupling." }));

    REGISTER(coupling_dummy, ECFMetaData()
            .setdescription({ "Settings for coupling dummy." }));

    dryrun = 0;
    REGISTER(dryrun, ECFMetaData()
            .setdescription({ "Just for testing assembler" })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));
    simple_output = false;
    REGISTER(simple_output, ECFMetaData()
            .setdescription({ "Compute solution independent on conductivity and elasticity." })
            .setdatatype({ ECFDataType::BOOL }));
    always_update_conductivity = false;
    REGISTER(always_update_conductivity, ECFMetaData()
            .setdescription({ "Always call update conductivity function." })
            .setdatatype({ ECFDataType::BOOL }));

    physics = PhysicsConfiguration::TYPE::HEAT_TRANSFER;
    REGISTER(physics, ECFMetaData()
            .setdescription({ "Physics" })
            .setdatatype({ ECFDataType::OPTION })
//            .addoption(ECFOption().setname("THERMO_ELASTICITY").setdescription("Coupled physics."))
            .addoption(ECFOption().setname("HEAT_TRANSFER").setdescription("Heat transfer."))
            .addoption(ECFOption().setname("STRUCTURAL_MECHANICS").setdescription("Structural mechanics 2D."))
//            .addoption(ECFOption().setname("ACOUSTICS").setdescription("Acoustics."))
            );

//    REGISTER(thermo_elasticity, ECFMetaData()
//            .setdescription({ "Coupled physics" })
//            .allowonly([&] () { return physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY; }));
    REGISTER(heat_transfer, ECFMetaData()
            .setdescription({ "Heat transfer" })
            .allowonly([&] () { return physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER; }));
    REGISTER(structural_mechanics, ECFMetaData()
            .setdescription({ "Structural mechanics" })
            .allowonly([&] () { return physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS; }));
//    REGISTER(acoustics, ECFMetaData()
//            .setdescription({ "Acoustics" })
//            .allowonly([&] () { return physics == PhysicsConfiguration::TYPE::ACOUSTICS; }));

    REGISTER(output, ECFMetaData()
            .setdescription({ "Output configurations." }));
}

ECF::ECF()
: mesh_morphing(this),
  output(this)
{
    info::ecf = this;
    _init();
}

ECF::ECF(const std::string &file)
: ECF()
{
    fill(file);
}

ECF::ECF(int *argc, char ***argv, const std::string &app)
: ECF()
{
    fill(argc, argv, app);
}

void ECF::fill(const std::string &file)
{
    ecffile = file;
    ECFReader::read(*this->ecfdescription, file, this->ranges, this->default_args, this->variables);
    set();
}

void ECF::fill(int *argc, char ***argv, const std::string &app)
{
    ecffile = ECFReader::read(*this->ecfdescription, argc, argv, this->ranges, this->default_args, this->variables);
    exe = std::string(info::system::buildpath()) + "/" + app;
    if (StringCompare::caseInsensitiveEq("dummycoupler", app)) {
        info::ecf->output.path = "dummy";
    }
    set();
}

void ECF::set()
{
    size_t namebegin = ecffile.find_last_of("/") + 1;
    size_t nameend = ecffile.find_last_of(".");
    name = ecffile.substr(namebegin, nameend - namebegin);

    struct tm *timeinfo;
    char datetime[80];
    timeinfo = std::localtime(&TimeLogger::initTime);
    std::strftime(datetime, 80, "%F-at-%Hh-%Mm-%Ss", timeinfo);

    outpath = info::ecf->output.path + "/" + std::string(datetime);

    if (info::mpi::grank) {
        Communication::barrier(MPITools::global);
    } else {
        utils::createDirectory(info::ecf->outpath + "/PREPOSTDATA");
        std::string symlink = info::ecf->output.path + "/last";
        if (utils::exists(symlink)) {
            utils::remove(symlink);
        }
        utils::createSymlink(datetime, symlink);
        utils::copyFile(info::ecf->ecffile, symlink + "/" + info::ecf->name + ".ecf");
        Communication::barrier(MPITools::global);
    }
}



