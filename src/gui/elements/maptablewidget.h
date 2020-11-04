
#ifndef MAPTABLEWIDGET_H
#define MAPTABLEWIDGET_H

#include "isavableobject.h"

#include "gui/declarations/tablewidget.h"

namespace espreso
{

class MapTableWidget : public TableWidget, public ISavableObject
{
    Q_OBJECT

public:
    virtual QVector<QPair<QString, QString> > dataInRows() = 0;

    using TableWidget::addData;
    virtual void addData(const QString& data) override;
    virtual QString data() override;

protected:
    MapTableWidget(const QStringList& headlines,
                   QWidget *parent = 0,
                   std::unique_ptr<CleanRowFactory> rowFactory = std::unique_ptr<CleanRowFactory>(new DefaultCleanRowFactory));
};

}

#endif // MAPTABLEWIDGET_H
