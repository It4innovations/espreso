#include "defaultmaptablewidget.h"

#include "validators/validatordelegate.h"
#include "validators/validatorfactory.h"

using namespace espreso;

DefaultMapTableWidget::DefaultMapTableWidget(ECFObject* map,
                                           const QStringList& headlines,
                                           QWidget* parent,
                                           std::unique_ptr<CleanRowFactory> rowFactory) :
    MapTableWidget(headlines, parent, std::move(rowFactory))
{
    this->m_map = map;

    ECFDataType fst = map->metadata.datatype[0];
    ECFDataType snd = map->metadata.datatype[1];
    this->setupColumn(0, fst);
    this->setupColumn(1, snd);
}

void DefaultMapTableWidget::setupColumn(int col, ECFDataType type)
{
    if (type == ECFDataType::NONNEGATIVE_INTEGER)
    {
        this->mTable->setItemDelegateForColumn(
                    col,
                    new ValidatorDelegate(new NonnegativeIntegerValidatorFactory(), this));
    }
    else if (type == ECFDataType::POSITIVE_INTEGER)
    {
        this->mTable->setItemDelegateForColumn(
                    col,
                    new ValidatorDelegate(new PositiveIntegerValidatorFactory(), this));
    }
}

QVector<QPair<QString, QString> > DefaultMapTableWidget::dataInRows()
{
    QVector<QPair<QString, QString> > ret;

    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex key_index = this->mModel->index(row, 0);
        QModelIndex val_index = this->mModel->index(row, 1);

        ret.append(QPair<QString, QString>(key_index.data().toString(),
                                           val_index.data().toString()));
    }

    return ret;
}

bool DefaultMapTableWidget::isValid()
{
    return true;
}

QString DefaultMapTableWidget::errorMessage()
{
    return QLatin1String("DefaultMapTableWidget: Error.");
}

void DefaultMapTableWidget::save()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex key_index = this->mModel->index(row, 0);
        QModelIndex val_index = this->mModel->index(row, 1);

        this->m_map->getParameter(key_index.data().toString().toStdString())
                ->setValue(val_index.data().toString().toStdString());
    }
}
