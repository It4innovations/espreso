#ifndef BOOLMAPTABLEWIDGET_H
#define BOOLMAPTABLEWIDGET_H

#include "defaultmaptablewidget.h"

#include "config/configuration.h"

namespace espreso
{

class BoolCleanRowFactory : public CleanRowFactory
{
public:
    QList<QStandardItem*> create(int columns) override;
};

class BoolMapTableWidget : public DefaultMapTableWidget
{
    Q_OBJECT
public:
    BoolMapTableWidget(ECFObject* map, const QStringList& headlines, QWidget* parent = 0);

    void addRow(const QVector<QString> &rowData) override;
    QVector<QPair<QString, QString> > dataInRows() override;
    QString errorMessage() override;
    void save() override;

private:
    QString checkBoxValue(int row);
};

}

#endif // BOOLMAPTABLEWIDGET_H
