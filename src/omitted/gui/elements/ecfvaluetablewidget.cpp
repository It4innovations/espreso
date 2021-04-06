#include "ecfvaluetablewidget.h"
#include "ui_ecfvaluetablewidget.h"

#include <QLineEdit>
#include "validators/validatorfactory.h"
#include "validators/validatordelegate.h"
#include "textitemwidgetfactory.h"
#include "textitemdelegate.h"

#include <QDebug>

using namespace espreso;

ECFValueTableWidget::ECFValueTableWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ECFValueTableWidget)
{
    ui->setupUi(this);

    this->m_model = new QStandardItemModel(ui->table);
    QStringList header;
    header << tr("Property") << tr("Value");
    this->m_model->setHorizontalHeaderLabels(header);
    ui->table->setModel(this->m_model);
    ui->table->setEditTriggers(QAbstractItemView::AllEditTriggers);

    connect(this->m_model, SIGNAL(itemChanged(QStandardItem*)),
            this, SLOT(onItemChanged(QStandardItem*)));
}


ECFValueTableWidget::~ECFValueTableWidget()
{
    delete ui;
}

void ECFValueTableWidget::add(ECFValue *value)
{
    if (value->metadata.datatype.empty())
    {
        qCritical("ECFValueTableWidget: Attempt to add an ECFValue with no data type!");
        return;
    }

    switch(value->metadata.datatype.at(0))
    {
    case ECFDataType::BOUNDARY_REGION:
    case ECFDataType::STRING:
        this->addString(value);
        break;
    case ECFDataType::NONNEGATIVE_INTEGER:
        this->addNonnegativeInteger(value);
        break;
    case ECFDataType::POSITIVE_INTEGER:
        this->addPositiveInteger(value);
        break;
    case ECFDataType::FLOAT:
        this->addFloat(value);
        break;
    case ECFDataType::ENUM_FLAGS:
    case ECFDataType::OPTION:
        this->addOption(value);
        break;
    case ECFDataType::BOOL:
        this->addBool(value);
        break;
    default: ;
    }

}

void ECFValueTableWidget::addString(ECFValue *string)
{
    ECFValue* value = string;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);
}

void ECFValueTableWidget::addNonnegativeInteger(ECFValue *nnint)
{
    ECFValue* value = nnint;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);

    ui->table->setItemDelegateForRow(this->m_model->rowCount() - 1,
                                     new ValidatorDelegate(new NonnegativeIntegerValidatorFactory));
}

void ECFValueTableWidget::addPositiveInteger(ECFValue *pint)
{
    ECFValue* value = pint;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);

    ui->table->setItemDelegateForRow(this->m_model->rowCount() - 1,
                                     new ValidatorDelegate(new PositiveIntegerValidatorFactory));
}

void ECFValueTableWidget::addFloat(ECFValue *p_float)
{
    ECFValue* value = p_float;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);

    ui->table->setItemDelegateForRow(this->m_model->rowCount() - 1,
                                     new ValidatorDelegate(new DoubleValidatorFactory));
}

void ECFValueTableWidget::addOption(ECFValue *option)
{
    ECFValue* value = option;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QString text = QString::fromStdString(value->getValue());
    QStandardItem* val = new QStandardItem(text);
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);

    QStringList options;
    for (auto it = value->metadata.options.cbegin();
         it != value->metadata.options.cend();
         ++it)
    {
        options << QString::fromStdString(it->name);
    }
    ui->table->setItemDelegateForRow(this->m_model->rowCount() - 1,
                                     new TextItemDelegate(text, new OptionWidgetFactory(options)));

}

void ECFValueTableWidget::addBool(ECFValue *p_bool)
{
    ECFValue* value = p_bool;
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);

    QStandardItem* val = new QStandardItem("");
    val->setEditable(false);
    val->setCheckable(true);
    if (value->getValue().compare("TRUE") == 0)
        val->setCheckState(Qt::Checked);
    else
        val->setCheckState(Qt::Unchecked);

    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);
}

void ECFValueTableWidget::addWithDelegate(ECFValue* value, QItemDelegate* delegate)
{
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);
    ui->table->setItemDelegateForRow(this->m_model->rowCount() - 1,
                                     delegate);
}

void ECFValueTableWidget::addWithWidget(ECFValue* value, QWidget *widget)
{
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    list << name << val;
    this->m_model->appendRow(list);
    ui->table->setIndexWidget(this->m_model->index(m_model->rowCount() - 1, 1),
                              widget);
}

void ECFValueTableWidget::resizeCellsToContent()
{
    ui->table->resizeColumnsToContents();
    ui->table->resizeRowsToContents();

    ui->table->horizontalHeader()->setStretchLastSection(true);
}

void ECFValueTableWidget::onItemChanged(QStandardItem *item)
{
    if (item->isCheckable()
            || this->m_values[item->row()]->metadata.datatype.at(0) == ECFDataType::OPTION)
    {
        emit itemChanged();
    }
}

bool ECFValueTableWidget::isValid()
{
    for (int i = 0; i < this->m_model->rowCount(); i++)
    {
        if (this->m_model->item(i, 1)->text().isEmpty()
                && !this->m_model->item(i, 1)->isCheckable())
        {
            this->m_invalid = i;
            return false;
        }
    }

    this->m_invalid = -1;
    return true;
}

QString ECFValueTableWidget::errorMessage()
{
    if (this->m_invalid != -1)
        return tr("Empty property %1").arg(this->m_model->item(m_invalid, 0)->text());
    return QLatin1String("");
}

void ECFValueTableWidget::save()
{
    for (int i = 0; i < this->m_values.size(); i++)
    {
        QStandardItem* item = this->m_model->item(i, 1);
        if (item->isCheckable())
        {
            if (item->checkState() == Qt::Checked)
                this->m_values[i]->setValue("TRUE");
            else
                this->m_values[i]->setValue("FALSE");
        }
        else
        {
            this->m_values[i]->setValue(
                        this->m_model->item(i, 1)->text().toStdString()
                        );
        }
    }
}

void ECFValueTableWidget::saveState()
{
    this->m_stored_values.clear();

    for (int i = 0; i < this->m_values.size(); i++)
    {
        this->m_stored_values << this->m_values[i]->getValue();
    }
}

void ECFValueTableWidget::restoreState()
{
    for (int i = 0; i < this->m_stored_values.size(); i++)
    {
        this->m_values[i]->setValue(this->m_stored_values[i]);
        QStandardItem *item = this->m_model->item(i, 1);
        if (item->isCheckable())
        {
            if (this->m_values[i]->getValue().compare("TRUE") == 0)
                item->setCheckState(Qt::Checked);
            else
                item->setCheckState(Qt::Unchecked);
        }
        else
        {
            item->setText(QString::fromStdString(this->m_values[i]->getValue()));
        }
    }
}
