
#ifndef SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_

#include "general/subkernel.h"
#include "basis/evaluator/evaluator.h"
#include "esinfo/stepinfo.h"
#include "math/primitives/vector_sparse.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct ECFExpression;
struct ECFExpressionVector;
struct ECFExpressionOptionalVector;
struct ECFHarmonicExpressionVector;
struct ConvectionConfiguration;
struct ImpedanceConfiguration;
struct PointSourceConfiguration;
struct PhysicsConfiguration;

class Assembler
{
public:
    Assembler(PhysicsConfiguration &settings);
    virtual ~Assembler();

    static ECFExpression*               getExpression(size_t interval        , std::map<std::string, ECFExpression> &settings);
    static ECFExpressionVector*         getExpression(size_t interval        , std::map<std::string, ECFExpressionVector> &settings);
    static ECFHarmonicExpressionVector* getExpression(size_t interval        , std::map<std::string, ECFHarmonicExpressionVector> &settings);
    static ECFExpression*               getExpression(const std::string &name, std::map<std::string, ECFExpression> &settings);
    static ECFExpressionVector*         getExpression(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);
    static ECFHarmonicExpressionVector* getExpression(const std::string &name, std::map<std::string, ECFHarmonicExpressionVector> &settings);

    static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpression> &settings);
    static Evaluator* getEvaluator(size_t interval, std::map<std::string, ECFExpressionVector> &settings, int dim);

    PhysicsConfiguration &settings;

protected:
    void assemble(const step::Step &step, SubKernel::Action action);
    virtual void elements(const step::Step &step, SubKernel::Action action, size_t interval) =0;
    virtual void boundary(const step::Step &step, SubKernel::Action action, size_t region, size_t interval) =0;
    virtual void nodes(const step::Step &step, SubKernel::Action action, size_t region, size_t interval) =0;
    virtual void bem(const step::Step &step, SubKernel::Action action, size_t domain, double *BETI)
    {
        eslog::error("implement BEM assembler\n");
    }

    bool checkExpression(const std::string &name, ECFExpression &expression);
    bool checkExpression(const std::string &name, ECFExpressionVector &expression);
    bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings);
    bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);
    bool checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, int dim);
    bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings);
    bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings);
    bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings);
    bool checkBoundaryParameter(const std::string &name, std::map<std::string, ECFHarmonicExpressionVector> &settings);

    void printElementVolume(std::vector<double> &volume);
    void printBoundarySurface(std::vector<double> &surface);

    void printMaterials(const std::map<std::string, std::string> &settings);
    template<typename Ttype>
    void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

    bool isBEM(size_t interval);

    std::vector<int> BEM;
    std::vector<double*> BETI;
    bool threaded;
};

template <typename T>
static bool reset(T *t, bool constant)
{
    if (t) {
        if (!t->filled || !constant) {
            t->set(0);
            t->updated = true;
            return true;
        }
        t->updated = false;
        return false;
    }
    return false;
}

template <typename T>
static void update(T *t, bool constant)
{
    if (t && t->updated) {
        t->filled = true;
        t->synchronize();
    }
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_ */
