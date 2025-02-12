
#include "transformation.h"
#include "config/configuration.hpp"

using namespace espreso;

InputTransformationConfiguration::InputTransformationConfiguration()
{
    transformation = TRANSFORMATION::TRANSLATE;
    REGISTER(transformation, ECFMetaData()
            .setdescription({ "A type of the transformation." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("TRANSLATE").setdescription("Move all nodes with defined direction."))
            .addoption(ECFOption().setname("ROTATE").setdescription("Rotate all nodes around ."))
            .addoption(ECFOption().setname("SCALE").setdescription("Use POSIX API."))
            .addoption(ECFOption().setname("SHEAR").setdescription("Use MPI_File_read_at."))
            );

    x = 0;
    REGISTER(x, ECFMetaData()
            .setdescription({ "X" })
            .setdatatype({ ECFDataType::FLOAT }));
    y = 0;
    REGISTER(y, ECFMetaData()
            .setdescription({ "Y" })
            .setdatatype({ ECFDataType::FLOAT }));
    z = 0;
    REGISTER(z, ECFMetaData()
            .setdescription({ "Z" })
            .setdatatype({ ECFDataType::FLOAT }));

    instances = 1;
    REGISTER(instances, ECFMetaData()
            .setdescription({ "N" })
            .setdatatype({ ECFDataType::INTEGER }));
}
