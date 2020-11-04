
#ifndef UTILS_H
#define UTILS_H

#include <QVector>

#include <string>

namespace espreso
{

class Utils
{
public:
    Utils();

    static int caseInsensitiveIndexOf(const QVector<std::string>& where, const std::string& what);
};

}

#endif // UTILS_H
