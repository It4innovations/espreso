
#ifndef SRC_ESINFO_MPIINFO_H_
#define SRC_ESINFO_MPIINFO_H_

#include "mpi.h"

namespace espreso {
namespace info {
namespace mpi {
    // instance MPI communication settings
    extern int rank;
    extern int size;
    extern MPI_Comm comm;

    // inter-instances MPI communication settings
    extern int irank;
    extern int isize;
    extern MPI_Comm icomm;

    // global MPI communication settings
    extern int grank;
    extern int gsize;
    extern MPI_Comm gcomm;

    // MPI_Init_thread provided level
    extern int threading;

    void init(int *argc, char ***argv);
    void init(MPI_Comm comm);
    bool divide(int meshDuplication);
    void finish();

    void print();
}
}
}



#endif /* SRC_ESINFO_MPIINFO_H_ */
