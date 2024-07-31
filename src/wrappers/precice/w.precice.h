
#ifndef SRC_WRAPPERS_PRECICE_W_PRECICE_H_
#define SRC_WRAPPERS_PRECICE_W_PRECICE_H_

namespace espreso {

struct PreciceData;

struct Precice {

    static void dummy();

    Precice(const char *name, bool active);
    ~Precice();

    void read(const char *name, double *data, double dt);
    void write(const char *name, double *data);
    void advance(double dt);

private:
    PreciceData *_data;
};

}

#endif /* SRC_WRAPPERS_PRECICE_W_PRECICE_H_ */
