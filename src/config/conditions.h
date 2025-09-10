
#ifndef SRC_CONFIG_CONDITIONS_H_
#define SRC_CONFIG_CONDITIONS_H_

#include <string>

namespace espreso {

struct ECFParameter;

struct GeneralValue {
    virtual bool equal(const void *data) const { return true; }
    virtual std::string tostring() const { return "TRUE"; }

    virtual ~GeneralValue() {}
    virtual GeneralValue* copy() const { return new GeneralValue(); }
};

struct EnumValue: public GeneralValue {
    virtual int index() const =0;
};

template <typename TValue>
struct EnumValueHolder: public EnumValue {
    TValue value;

    bool equal(const void *data) const { return *static_cast<const TValue*>(data) == value; }
    int index() const { return static_cast<int>(value); }

    EnumValueHolder(const TValue &value): value(value) {}
    GeneralValue* copy() const { return new EnumValueHolder<TValue>(value); }
};

template <typename TValue>
struct GeneralValueHolder: public GeneralValue {
    TValue value;

    bool equal(const void *data) const { return *static_cast<const TValue*>(data) == value; }
    std::string tostring() const { return std::to_string(value); }

    GeneralValueHolder(const TValue &value): value(value) {}
    GeneralValue* copy() const { return new GeneralValueHolder<TValue>(value); }
};

template<>
inline std::string GeneralValueHolder<bool>::tostring() const
{
    if (this->value) return "1";
    else return "0";
}


struct ECFAbstractCondition
{

public:
    virtual bool evaluate() const = 0;

    virtual bool isset() const = 0;
    virtual bool match(const void* parameter) const = 0;
    virtual void bind(
        const ECFParameter* parameter, 
        const std::string& prefix = ""
    ) = 0;
    virtual std::string compose() const = 0;

    virtual ECFAbstractCondition* copy() const = 0;

    virtual ~ECFAbstractCondition() {}
};

struct ECFFalseCondition : public ECFAbstractCondition
{

public:
    bool evaluate() const override { return false; }

    bool isset() const override { return true; }
    bool match(const void*) const override { return true; }
    void bind(
        const ECFParameter*,
        const std::string& = ""
    ) override {}
    std::string compose() const override { return "false"; }

    ECFAbstractCondition* copy() const override { return new ECFFalseCondition(); }
};

struct ECFCondition : public ECFAbstractCondition
{
protected:
    const void *parameter;
    const int operation;
    const GeneralValue *value;
    std::string prefix;
    bool bound;
    const ECFParameter* ecfparameter;

public:
    static const int EQUALS = 0;
    static const int NOT_EQUALS = 1;


    bool evaluate() const override;

    bool isset() const override { return parameter; }
    bool match(const void* parameter) const override 
        { return this->parameter == parameter; }
    void bind(
        const ECFParameter* parameter, 
        const std::string& prefix = ""
    ) override;
    std::string compose() const override;

    ECFCondition()
    : parameter(NULL), operation(EQUALS), value(new GeneralValue()),
    prefix(""), bound(true), ecfparameter(NULL) {}

    template <typename TValue>
    typename std::enable_if<std::is_enum<TValue>::value, const GeneralValue*>::type
    init(const TValue &value) { return new EnumValueHolder<TValue>(value); }

    template <typename TValue>
    typename std::enable_if<!std::is_enum<TValue>::value, const GeneralValue*>::type
    init(const TValue &value) { return new GeneralValueHolder<TValue>(value); }

    template <typename TParameter, typename TValue>
    ECFCondition(const TParameter &parameter, const int operation, const TValue &value)
    : parameter(&parameter), operation(operation), value(init<TValue>(value)), bound(false),
    ecfparameter(NULL) {}

    ECFCondition(const ECFCondition &other)
    : parameter(other.parameter), operation(other.operation), 
    value(other.value->copy()), prefix(other.prefix), bound(other.bound) 
    { if (other.bound) this->ecfparameter = other.ecfparameter; }

    ECFCondition(ECFCondition &&other)
    : parameter(std::move(other.parameter)), operation(other.operation), 
    value(other.value->copy()), prefix(other.prefix), bound(other.bound)
    { if (other.bound) this->ecfparameter = other.ecfparameter; }

    ECFCondition* copy() const override { return new ECFCondition(*this); }

    virtual ~ECFCondition() { delete value; }
};

struct ECFConditionPair : public ECFAbstractCondition
{
protected:
    ECFAbstractCondition* left;
    const int operation;
    ECFAbstractCondition* right;

public:
    static const int AND = 0;
    static const int OR = 1;

    bool evaluate() const override;

    bool isset() const override { return true; }
    bool match(const void* parameter) const override
        {return this->left->match(parameter) || this->right->match(parameter);}
    void bind(
        const ECFParameter* parameter, 
        const std::string& prefix = ""
    ) override;
    std::string compose() const override;

    ECFConditionPair(const ECFAbstractCondition& leftOperand,
        const int operation,
        const ECFAbstractCondition& rightOperand) 
        : left(leftOperand.copy()), operation(operation), right(rightOperand.copy()) {}

    ECFConditionPair(const ECFConditionPair &other)
    : left(other.left->copy()), operation(other.operation), right(other.right->copy()) {}

    ECFConditionPair(ECFConditionPair &&other)
    : left(std::move(other.left)), operation(other.operation), right(std::move(other.right)) 
    {
        other.left = NULL;
        other.right = NULL;
    }
    
    ECFAbstractCondition* copy() const override { return new ECFConditionPair(*this); };

    virtual ~ECFConditionPair();
};

inline ECFConditionPair operator&(const ECFAbstractCondition& c1, const ECFAbstractCondition& c2)
{
    return {c1, ECFConditionPair::AND, c2};
}

inline ECFConditionPair operator|(const ECFAbstractCondition& c1, const ECFAbstractCondition& c2)
{
    return {c1, ECFConditionPair::OR, c2};
}

}

#endif /* SRC_CONFIG_CONDITIONS_H_ */
