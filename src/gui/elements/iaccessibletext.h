
#ifndef IACCESSIBLETEXT_H
#define IACCESSIBLETEXT_H

#include <QString>

namespace espreso
{

class IAccessibleText
{
public:
    virtual ~IAccessibleText() {}

    virtual void setText(const QString& text) = 0;
    virtual QString text() = 0;
};

}

#endif // IACCESSIBLETEXT_H
