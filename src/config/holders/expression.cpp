
#include "expression.h"
#include "config/configuration.hpp"
#include "esinfo/eslog.hpp"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "basis/evaluator/constevaluator.h"
#include "basis/evaluator/expressionevaluator.h"
#include "basis/evaluator/tableinterpolationevaluator.h"

#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

bool ECFExpression::forall(const std::map<std::string, ECFExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
    for (auto it = parameter.begin(); it != parameter.end(); ++it) {
        if (!fnc(it->second)) {
            return false;
        }
    }
    return true;
}

bool ECFHarmonicExpression::forall(const std::map<std::string, ECFHarmonicExpression> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
    for (auto it = parameter.begin(); it != parameter.end(); ++it) {
        if (!fnc(it->second.magnitude)) {
            return false;
        }
        if (!fnc(it->second.phase)) {
            return false;
        }
    }
    return true;
}

bool ECFExpressionVector::forall(const std::map<std::string, ECFExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
    for (auto it = parameter.begin(); it != parameter.end(); ++it) {
        if (!fnc(it->second.x) || !fnc(it->second.y) || !fnc(it->second.z)) {
            return false;
        }
    }
    return true;
}

bool ECFHarmonicExpressionVector::forall(const std::map<std::string, ECFHarmonicExpressionVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
    for (auto it = parameter.begin(); it != parameter.end(); ++it) {
        if (!fnc(it->second.magnitude.x) || !fnc(it->second.magnitude.y) || !fnc(it->second.magnitude.z)) {
            return false;
        }
        if (!fnc(it->second.phase.x) || !fnc(it->second.phase.y) || !fnc(it->second.phase.z)) {
            return false;
        }
    }
    return true;
}

bool ECFExpressionOptionalVector::forall(const std::map<std::string, ECFExpressionOptionalVector> &parameter, std::function<bool(const ECFExpression &expr)> fnc)
{
    for (auto it = parameter.begin(); it != parameter.end(); ++it) {
        if (it->second.all.isset) {
            if (!fnc(it->second.all)) {
                return false;
            }
        } else {
            if (!fnc(it->second.x) || !fnc(it->second.y) || !fnc(it->second.z)) {
                return false;
            }
        }
    }
    return true;
}

ECFExpression::ECFExpression()
: evaluator(NULL), isset(false)
{
    evaluator = Evaluator::create(value);
}

ECFExpression::~ECFExpression()
{
    if (evaluator) {
        delete evaluator;
    }
}

ECFHarmonicExpression::ECFHarmonicExpression()
: type(Type::COMPONENTS)
{
    type = Type::COMPONENTS;
    REGISTER(type, ECFMetaData()
            .setdescription({ "Description type." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("COMPONENTS").setdescription("Separated components.")));

    REGISTER(magnitude, ECFMetaData().setdescription({ "Magnitude." }));
    REGISTER(phase, ECFMetaData().setdescription({ "Phase." }));
}

ECFHarmonicExpressionVector::ECFHarmonicExpressionVector()
: type(Type::COMPONENTS)
{
    type = Type::COMPONENTS;
    REGISTER(type, ECFMetaData()
            .setdescription({ "Description type." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("COMPONENTS").setdescription("Separated components.")));

    REGISTER(magnitude, ECFMetaData().setdescription({ "Magnitude." }));
    REGISTER(phase, ECFMetaData().setdescription({ "Phase." }));
}

ECFExpressionVector::ECFExpressionVector()
{
    REGISTER(x, ECFMetaData()
            .setdescription({ "x-direction." })
            .setdatatype({ ECFDataType::EXPRESSION }));
    REGISTER(y, ECFMetaData()
            .setdescription({ "y-direction." })
            .setdatatype({ ECFDataType::EXPRESSION }));
    REGISTER(z, ECFMetaData()
            .setdescription({ "z-direction." })
            .setdatatype({ ECFDataType::EXPRESSION }));
}

ECFExpressionOptionalVector::ECFExpressionOptionalVector()
{
    REGISTER(all, ECFMetaData()
            .setdescription({ "all-directions." })
            .setdatatype({ ECFDataType::EXPRESSION }));
}


