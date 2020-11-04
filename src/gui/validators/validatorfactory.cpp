
#include "validatorfactory.h"

QValidator* DoubleValidatorFactory::create(QObject *parent) const
{
    return new QRegExpValidator(QRegExp("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$"), parent);
}

QValidator* PositiveIntegerValidatorFactory::create(QObject *parent) const
{
    return new QRegExpValidator(QRegExp("^[1-9][0-9]*$"), parent);
}

QValidator* NonnegativeIntegerValidatorFactory::create(QObject *parent) const
{
    return new QRegExpValidator(QRegExp("^[0]?([1-9][0-9]*)$"), parent);
}
