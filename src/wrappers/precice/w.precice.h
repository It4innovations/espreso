
#ifndef SRC_WRAPPERS_PRECICE_W_PRECICE_H_
#define SRC_WRAPPERS_PRECICE_W_PRECICE_H_

namespace espreso {

struct PreciceData;

struct Precice {

    static void dummy();

    Precice();
    ~Precice();

    double timeStep(double dt);

    bool requiresWritingCheckpoint();
    bool requiresReadingCheckpoint();

    void read(double *data, double dt);
    void write(double *data);
    void advance(double dt);

private:
    PreciceData *_data;
};

}

#endif /* SRC_WRAPPERS_PRECICE_W_PRECICE_H_ */
