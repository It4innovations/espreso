
#ifndef SRC_BASIS_EVALUATOR_EVALUATOR_H_
#define SRC_BASIS_EVALUATOR_EVALUATOR_H_

#include "parameter.h"

#include <vector>

namespace espreso {

class Evaluator {
public:
    static Evaluator* create(const std::string &expression);

    Evaluator();
    virtual ~Evaluator() {};
    virtual double evaluate(int t = 0) const { return 0; }

    double& getParameter(const std::string &name, int t = 0);

    double& getSubstep(int t = 0)     { return getParameter("SUBSTEP", t); }
    double& getTime(int t = 0)        { return getParameter("TIME", t); }
    double& getFrequency(int t = 0)   { return getParameter("FREQUENCY", t); }

    double& getX(int t = 0) { return getParameter("X", t); }
    double& getY(int t = 0) { return getParameter("Y", t); }
    double& getZ(int t = 0) { return getParameter("Z", t); }
    double& getCoordinateX(int t = 0) { return getParameter("COORDINATE_X", t); }
    double& getCoordinateY(int t = 0) { return getParameter("COORDINATE_Y", t); }
    double& getCoordinateZ(int t = 0) { return getParameter("COORDINATE_Z", t); }

    double& getTemperature(int t = 0) { return getParameter("TEMPERATURE", t); }

    bool isConst() const { return parameters.front().size() == 1; }
    bool needTemperature(int t = 0) { return &getTemperature(t) != &parameters[t].front().value; }
    bool needCoordinates(int t = 0) { return &getCoordinateX(t) != &parameters[t].front().value || &getCoordinateY(t) != &parameters[t].front().value || &getCoordinateZ(t) != &parameters[t].front().value; }

    std::vector<std::vector<EvaluatorParameter> > parameters;
    virtual const char* expression() const { return "null"; }
};

}



#endif /* SRC_BASIS_EVALUATOR_EVALUATOR_H_ */
