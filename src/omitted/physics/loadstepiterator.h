
#ifndef SRC_PHYSICS_LOADSTEPITERATOR_H_
#define SRC_PHYSICS_LOADSTEPITERATOR_H_

namespace espreso {

class LoadStepSolver;
class SubStepSolver;
class LinearSystem;

class LoadStepIterator {

public:
    LoadStepIterator();
    ~LoadStepIterator();

    void prepareExpressions();
    bool next();

protected:

    template <typename TPhysics> bool next(TPhysics &configuration);

    LoadStepSolver *_loadStepSolver;
    SubStepSolver *_subStepSolver;
    LinearSystem *_system;
};

}

#endif /* SRC_PHYSICS_LOADSTEPITERATOR_H_ */
