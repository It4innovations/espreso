
#include "utils.h"

using namespace espreso;

Utils::Utils()
{

}

int Utils::caseInsensitiveIndexOf(const QVector<std::string>& where, const std::string& what)
{
    int index = -1;
    int current = 0;

    QString key = QString::fromStdString(what);

    for (auto it = where.cbegin(); it != where.cend(); ++it)
    {
        QString item = QString::fromStdString(*it);
        if (item.compare(key, Qt::CaseInsensitive) == 0)
        {
            return current;
        }
        current++;
    }

    return index;
}
