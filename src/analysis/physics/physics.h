
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; }

struct Physics {
    virtual ~Physics() {}

    virtual bool analyze(step::Step &step) =0;
    virtual bool run(step::Step &step, Physics *prev) =0;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
