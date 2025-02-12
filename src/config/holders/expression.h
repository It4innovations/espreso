
#ifndef SRC_CONFIG_HOLDERS_EXPRESSION_H_
#define SRC_CONFIG_HOLDERS_EXPRESSION_H_

#include "config/metadata.h"
#include "config/description.h"
#include <vector>
#include <map>
#include <functional>

namespace espreso {

class Evaluator;
struct BoundaryRegionStore;
struct ElementsRegionStore;

struct ECFExpression {
    std::string value;
    Evaluator *evaluator;
    bool isset;

    ECFExpression();
    ~ECFExpression();

    ECFExpression(const ECFExpression &other) = delete;
    ECFExpression& operator=(const ECFExpression &other) = delete;

    static bool forall(const std::map<std::string, ECFExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

struct ECFHarmonicExpression: public ECFDescription {
    enum class Type {
        COMPONENTS
    };

    Type type;
    ECFExpression magnitude, phase;

    ECFHarmonicExpression();

    static bool forall(const std::map<std::string, ECFHarmonicExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

struct ECFExpressionVector: public ECFDescription {
    ECFExpression data[3];
    ECFExpression &x = data[0], &y = data[1], &z = data[2];

    ECFExpressionVector(const ECFExpressionVector &other) = delete;
    ECFExpressionVector(const ECFExpressionVector &&other) = delete;
    ECFExpressionVector& operator=(const ECFExpressionVector &other) = delete;
    ECFExpressionVector& operator=(const ECFExpressionVector &&other) = delete;

    ECFExpressionVector();

    static bool forall(const std::map<std::string, ECFExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

struct ECFHarmonicExpressionVector: public ECFDescription {
    enum class Type {
        COMPONENTS
    };

    Type type;
    ECFExpressionVector magnitude, phase;

    ECFHarmonicExpressionVector();

    ECFHarmonicExpressionVector(const ECFHarmonicExpressionVector &other) = delete;
    ECFHarmonicExpressionVector& operator=(const ECFHarmonicExpressionVector &other) = delete;

    static bool forall(const std::map<std::string, ECFHarmonicExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

struct ECFExpressionOptionalVector: public ECFExpressionVector {
    ECFExpression all;

    ECFExpressionOptionalVector();

    ECFExpressionOptionalVector(const ECFExpressionOptionalVector &other) = delete;
    ECFExpressionOptionalVector& operator=(const ECFExpressionOptionalVector &other) = delete;

    static bool forall(const std::map<std::string, ECFExpressionOptionalVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc);
};

}

#endif /* SRC_CONFIG_HOLDERS_EXPRESSION_H_ */
