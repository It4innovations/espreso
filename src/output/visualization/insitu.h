
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_

#include "visualization.h"

namespace espreso {

class Catalyst;

struct InSitu: public Visualization {

    InSitu();
    ~InSitu();

    void updateMesh();
    void updateMonitors(step::TYPE type);
    void updateSolution(const step::Step &step, const step::Time &time);
    void updateSolution(const step::Step &step, const step::Frequency &frequency);

protected:
    Catalyst *_catalyst;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_ */
