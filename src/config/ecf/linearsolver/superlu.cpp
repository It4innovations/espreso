
#include "superlu.h"

#include <cmath>
#include "esinfo/mpiinfo.h"
#include "config/configuration.hpp"

using namespace espreso;

SuperLUConfiguration::SuperLUConfiguration()
{
    int procs_in_row = sqrt(info::mpi::size);

    np_row = procs_in_row;
    REGISTER(np_row, ECFMetaData()
        .setdescription({ "Number of processes per row of SuperLU grid of processes" })
        .setdatatype({ ECFDataType::INTEGER }));

    np_col = procs_in_row;
    REGISTER(np_col, ECFMetaData()
        .setdescription({ "Number of processes per column of SuperLU grid of processes" })
        .setdatatype({ ECFDataType::INTEGER }));
};
