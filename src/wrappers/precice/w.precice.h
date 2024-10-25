
#ifndef SRC_WRAPPERS_PRECICE_W_PRECICE_H_
#define SRC_WRAPPERS_PRECICE_W_PRECICE_H_

#include <string>

namespace espreso {

struct PreciceData;

struct Precice {

    Precice();
    ~Precice();

    double timeStep(double dt);

    bool requiresWritingCheckpoint();
    bool requiresReadingCheckpoint();

    void read(double dt);
    void write();
    void advance(double dt);

private:
    void _read(double *data, const std::string &name, double dt);
    void _write(double *data, const std::string &name);

    PreciceData *_data;
};

}

#endif /* SRC_WRAPPERS_PRECICE_W_PRECICE_H_ */
