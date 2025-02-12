
#ifndef SRC_WRAPPERS_PAPI_W_PAPI_H_
#define SRC_WRAPPERS_PAPI_W_PAPI_H_

namespace espreso {

class PAPI {

public:
    PAPI();
    ~PAPI();

    void init();
    void read(long long *values);

    static int status, values;
private:
    int set;
};

}

#endif /* SRC_WRAPPERS_PAPI_W_PAPI_H_ */
