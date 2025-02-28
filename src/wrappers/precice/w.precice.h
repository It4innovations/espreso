
#ifndef SRC_WRAPPERS_PRECICE_W_PRECICE_H_
#define SRC_WRAPPERS_PRECICE_W_PRECICE_H_

#include "config/ecf/coupling.h"

#include <string>

namespace espreso {

struct PreciceWrapper;

struct Precice {

    Precice(const CouplingConfiguration &configuration);
    ~Precice();

    double timeStep(double dt);

    bool requiresWritingCheckpoint();
    bool requiresReadingCheckpoint();

    void read(double dt);
    void write();
    void advance(double dt);

    void dummy();

private:
    void _read(double *data, const std::string &name, double dt);
    void _write(double *data, const std::string &name);

    PreciceWrapper *wrapper;
};

}

#endif /* SRC_WRAPPERS_PRECICE_W_PRECICE_H_ */
