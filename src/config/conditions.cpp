
#include "conditions.h"
#include "configuration.h"

#include <sstream>

using namespace espreso;

bool ECFCondition::evaluate() const
{
    if (!parameter) return true;

    switch(operation)
    {
    case NOT_EQUALS:
        return !value->equal(parameter);
    default:
        return value->equal(parameter);
    }
}

void ECFCondition::bind(
        const ECFParameter* parameter, 
        const std::string& prefix
)
{
    if (
        parameter->data() != this->parameter 
        || this->bound
    ) return;
    this->ecfparameter = parameter;
    this->prefix = prefix;
    this->bound = true;
}

std::string ECFCondition::compose() const
{
    if (!this->bound) return "";
    
    std::stringstream ss;

    // LEFT OPERAND
    if (!this->prefix.empty()) { ss << this->prefix << "."; }
    ss << ecfparameter->name;

    // OPERATOR
    switch (this->operation)
    {
    case EQUALS:
        ss << " == ";
        break;
    case NOT_EQUALS:
        ss << " != ";
        break;
    default:
        ss << " == ";
    }

    // RIGHT OPERAND
    switch (ecfparameter->metadata.datatype.front()) {
    case ECFDataType::ENUM_FLAGS:
        ss << "'"
        << ecfparameter->metadata
            .options[
                dynamic_cast<const EnumValue*>(value)->index() / 2
            ]
            .name
        << "'";
        break;
    case ECFDataType::OPTION:
        ss << "'"
        << ecfparameter->metadata
            .options[dynamic_cast<const EnumValue*>(value)->index()]
            .name
        << "'";
        break;
    default:
        ss << value->tostring();
    }

    return ss.str();
}

bool ECFConditionPair::evaluate() const
{
    bool leftEval = this->left->evaluate();
    bool rightEval = this->right->evaluate();

    if (this->operation == AND) return leftEval && rightEval;
    else return leftEval || rightEval;
}

void ECFConditionPair::bind(
        const ECFParameter* parameter, 
        const std::string& prefix
)
{
    this->left->bind(parameter, prefix);
    this->right->bind(parameter, prefix);
}

std::string ECFConditionPair::compose() const
{
    std::stringstream ss;
    // LEFT OPERAND
    ss << "(" << this->left->compose() << " ";
    // OPERATOR
    if (this->operation == AND) ss << "&&";
    else ss << "||";
    // RIGHT OPERAND
    ss << " " << this->right->compose() << ")";

    return ss.str();
}

ECFConditionPair::~ECFConditionPair()
{
    if (left) delete this->left; 
    if (right) delete this->right;
}
