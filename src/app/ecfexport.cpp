#include <iostream>
#include <fstream>
#include <sstream>

#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "config/export/json.h"
#include "config/ecf/ecf.h"

using namespace espreso;

int main(int argc, char **argv)
{
    info::system::setSignals();
    info::mpi::init(&argc, &argv);

    ECF ecf;
    if (info::mpi::rank == 0)
    {
        ECFJSONExport json(ecf.heat_transfer_2d.ecfdescription, std::cout);
        json.exportToStream();

        // ecf.forEachParameters([&] (const ECFParameter *parameter) {
        //     std::cout << parameter->name << std::endl;
        //     if (parameter->metadata.condition->isset()) {
        //         std::cout << parameter->name << ": ";
        //         ecf.forEachParameters([&] (const ECFParameter *conditionalParameter) {
        //             if (parameter->metadata.condition->match(conditionalParameter->data())) {
        //                 std::cout << parameter->metadata.condition->compose(conditionalParameter);
        //             }
        //         }, false, true);
        //         std::cout << "\n";
        //     }
        // }, false, true);
    } 

    info::mpi::finish();
    return 0;

    return 0;
}



