#include "boolmaptablewidget.h"

#include "validators/validatordelegate.h"
#include "validators/validatorfactory.h"

using namespace espreso;

QList<QStandardItem*> BoolCleanRowFactory::create(int columns)
{
    if (columns != 2) return QList<QStandardItem*>();

    QList<QStandardItem*> row;

    row << new QStandardItem();
    QStandardItem* check = new QStandardItem();
    check->setCheckable(true);
    check->setEditable(false);
    row << check;

    return row;
}

BoolMapTableWidget::BoolMapTableWidget(ECFObject* map,
                                           const QStringList& headlines,
                                           QWidget* parent) :
    DefaultMapTableWidget(map, headlines, parent, std::unique_ptr<CleanRowFactory>(new BoolCleanRowFactory))
{
}

void BoolMapTableWidget::addRow(const QVector<QString> &rowData)
{
    if (rowData.size() != 2)
    {
        qCritical("BoolMapTableWidget: Row should contain two values");
        return;
    }

    QList<QStandardItem*> row;
    row.append(new QStandardItem(rowData[0]));

    QStandardItem* check = new QStandardItem();
    check->setCheckable(true);
    check->setEditable(false);
    if (rowData[1].compare(QLatin1String("TRUE")) == 0)
        check->setCheckState(Qt::Checked);
    else check->setCheckState(Qt::Unchecked);

    this->mModel->insertRow(mModel->rowCount() - 1, row);
}

QString BoolMapTableWidget::checkBoxValue(int row)
{
    QStandardItem* val_index = this->mModel->item(row, 1);
    QString value;
    if (val_index->checkState() == Qt::Checked)
        value = QLatin1String("TRUE");
    else
        value = QLatin1String("FALSE");

    return value;
}

QVector<QPair<QString, QString> > BoolMapTableWidget::dataInRows()
{
    QVector<QPair<QString, QString> > ret;

    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QStandardItem* key_index = this->mModel->item(row, 0);
        QString value = this->checkBoxValue(row);

        ret.append(QPair<QString, QString>(key_index->data().toString(),
                                           value));
    }

    return ret;
}

QString BoolMapTableWidget::errorMessage()
{
    return QLatin1String("DefaultMapTableWidget: Error.");
}


void BoolMapTableWidget::save()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex key_index = this->mModel->index(row, 0);
        QString value = this->checkBoxValue(row);

        this->m_map->getParameter(key_index.data().toString().toStdString())
                ->setValue(value.toStdString());
    }
}
