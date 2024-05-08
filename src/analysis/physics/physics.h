
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; }

struct Physics {
    virtual ~Physics() {}

    virtual void analyze(step::Step &step) =0;
    virtual void run(step::Step &step) =0;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
