
#ifndef VALIDATORFACTORY_H
#define VALIDATORFACTORY_H

#include <QValidator>

class ValidatorFactory
{

public:
    ValidatorFactory() {}
    virtual ~ValidatorFactory() {}

    virtual QValidator* create(QObject* parent = nullptr) const = 0;
};


class DoubleValidatorFactory : public ValidatorFactory
{

public:
    QValidator* create(QObject *parent = nullptr) const override;
};


class PositiveIntegerValidatorFactory : public ValidatorFactory
{
public:
    QValidator* create(QObject *parent = nullptr) const override;
};

class NonnegativeIntegerValidatorFactory : public ValidatorFactory
{
public:
    QValidator* create(QObject *parent = nullptr) const override;
};

#endif // VALIDATORFACTORY_H
