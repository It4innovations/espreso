
#ifndef SRC_ANALYSIS_PATTERN_SYNCHRONIZATION_H_
#define SRC_ANALYSIS_PATTERN_SYNCHRONIZATION_H_

namespace espreso {

template <typename Object>
struct Synchronization {

    virtual ~Synchronization() {};

    virtual void gatherFromUpper(Object &o) =0;
    virtual void scatterToUpper(Object &o) =0;

    virtual void synchronize(Object &o)
    {
        gatherFromUpper(o);
        scatterToUpper(o);
    }
};

}

#endif /* SRC_ANALYSIS_PATTERN_SYNCHRONIZATION_H_ */
