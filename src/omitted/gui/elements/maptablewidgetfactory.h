#ifndef MAPTABLEWIDGETFACTORY_H
#define MAPTABLEWIDGETFACTORY_H

#include "maptablewidget.h"
#include "maptables/defaultmaptablewidget.h"
#include "maptables/boolmaptablewidget.h"

namespace espreso
{

class MapTableWidgetFactory
{
public:
    MapTableWidgetFactory();
    MapTableWidget* create(ECFObject* map, QWidget* parent = 0);

private:
    QVector<ECFDataType> m_first;
    QVector<ECFDataType> m_second;
};

}

#endif // MAPTABLEWIDGETFACTORY_H
