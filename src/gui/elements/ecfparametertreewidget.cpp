
#include "ecfparametertreewidget.h"
#include "ui_ecfparametertreewidget.h"
#include "textitemwidgetfactory.h"
#include "textitemdelegate.h"

#include "gui/validators/validatorfactory.h"
#include "gui/validators/validatordelegate.h"
#include "gui/declarations/datasetswidget.h"

#include <QLineEdit>
#include <QDebug>

using namespace espreso;

ECFParameterTreeWidget::ECFParameterTreeWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ECFParameterTreeWidget)
{
    ui->setupUi(this);

    this->m_model = new QStandardItemModel(ui->tree);
    this->m_root = this->m_model->invisibleRootItem();
    this->m_groups.append(this->m_root);

    QStringList header;
    header << tr("Property") << tr("Value");
    this->m_model->setHorizontalHeaderLabels(header);

    ui->tree->setModel(this->m_model);
    ui->tree->setEditTriggers(QAbstractItemView::AllEditTriggers);

    this->m_delegate = new ECFParameterTreeDelegate;
    ui->tree->setItemDelegateForColumn(1, this->m_delegate);

    connect(this->m_model, SIGNAL(itemChanged(QStandardItem*)),
            this, SLOT(onItemChanged(QStandardItem*)));
}

ECFParameterTreeWidget::~ECFParameterTreeWidget()
{
    delete ui;
}

int ECFParameterTreeWidget::createGroup(const QString& name, int parent_id)
{
    int id = this->m_groups.size();

    if (parent_id >= id)
    {
        qCritical("ECFParameterTreeWidget: Attempt to create a group with invalid ID of parent");
        return -1;
    }

    QList<QStandardItem*> row;
    QStandardItem* group = new QStandardItem(name);
    group->setEditable(false);
    QStandardItem* empty = new QStandardItem;
    empty->setEditable(false);
    row << group << empty;
    this->m_groups[parent_id]->appendRow(row);
    this->m_groups.append(group);
    if (parent_id == 0)
    {
        ui->tree->setExpanded(group->index(), true);
    }

    return id;
}

void ECFParameterTreeWidget::add(ECFValue *value, int parent_id)
{
    if (value->metadata.datatype.empty())
    {
        qCritical("ECFParameterTreeWidget: Attempt to add an ECFValue with no data type!");
        return;
    }

    switch(value->metadata.datatype.at(0))
    {
    case ECFDataType::REGION:
    case ECFDataType::STRING:
        this->addString(value, parent_id);
        break;
    case ECFDataType::NONNEGATIVE_INTEGER:
        this->addNonnegativeInteger(value, parent_id);
        break;
    case ECFDataType::POSITIVE_INTEGER:
        this->addPositiveInteger(value, parent_id);
        break;
    case ECFDataType::FLOAT:
        this->addFloat(value, parent_id);
        break;
    case ECFDataType::ENUM_FLAGS:
    case ECFDataType::OPTION:
        this->addOption(value, parent_id);
        break;
    case ECFDataType::BOOL:
        this->addBool(value, parent_id);
        break;
    default: ;
    }

}

void ECFParameterTreeWidget::addString(ECFValue *string, int parent_id)
{
    this->addRow(string,
                 new TextWidgetFactory,
                 parent_id);
}

void ECFParameterTreeWidget::addNonnegativeInteger(ECFValue *nnint, int parent_id)
{
    this->addRow(nnint,
                 new TextWidgetFactory(new NonnegativeIntegerValidatorFactory),
                 parent_id);
}

void ECFParameterTreeWidget::addPositiveInteger(ECFValue *pint, int parent_id)
{
    this->addRow(pint,
                 new TextWidgetFactory(new PositiveIntegerValidatorFactory),
                 parent_id);
}

void ECFParameterTreeWidget::addFloat(ECFValue *p_float, int parent_id)
{
    this->addRow(p_float,
                 new TextWidgetFactory(new DoubleValidatorFactory),
                 parent_id);
}

void ECFParameterTreeWidget::addOption(ECFValue *option, int parent_id)
{
    QStringList options;
    for (auto it = option->metadata.options.cbegin();
         it != option->metadata.options.cend();
         ++it)
    {
        options << QString::fromStdString(it->name);
    }

    this->addRow(option, new OptionWidgetFactory(options), parent_id);
}

void ECFParameterTreeWidget::addBool(ECFValue *p_bool, int parent_id)
{
    this->addRow(p_bool, nullptr, parent_id, true);
}

void ECFParameterTreeWidget::addWithEditorFactory(ECFValue* value,
                                                  TextItemWidgetFactory* editorFactory,
                                                  int parent_id)
{
    this->addRow(value, editorFactory, parent_id);
}

void ECFParameterTreeWidget::addRow(ECFValue* value, TextItemWidgetFactory* editorFactory,
                                    int parent_id, bool valueCheckable)
{
    int id = this->m_values.size();
    this->m_values.append(value);

    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);

    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    val->setData(id, ECFParameterTreeDelegate::EditorRole);

    if (valueCheckable)
    {
        val->setText("");
        val->setEditable(false);
        val->setCheckable(true);
        if (value->getValue().compare("TRUE") == 0)
            val->setCheckState(Qt::Checked);
        else
            val->setCheckState(Qt::Unchecked);
    }


    QList<QStandardItem*> list;
    list << name << val;

    this->m_groups[parent_id]->appendRow(list);
    this->m_delegate->registerEditor(id, editorFactory);
}

void ECFParameterTreeWidget::addWithWidget(ECFValue* value, QWidget *widget, int parent_id)
{
    int id = this->m_values.size();
    this->m_values.append(value);
    QStandardItem* name = new QStandardItem(
                QString::fromStdString(value->metadata.description.at(0))
                );
    name->setEditable(false);
    QStandardItem* val = new QStandardItem(
                QString::fromStdString(value->getValue())
                );
    QList<QStandardItem*> list;
    val->setData(id, ECFParameterTreeDelegate::EditorRole);
    list << name << val;
    this->m_groups[parent_id]->appendRow(list);
    ui->tree->setIndexWidget(val->index(),
                              widget);
}

void ECFParameterTreeWidget::resizeCellsToContent()
{
    ui->tree->resizeColumnToContents(0);
    ui->tree->resizeColumnToContents(1);
}

bool ECFParameterTreeWidget::findIndex(QStandardItem* item, QStandardItem *root, int *index)
{
    int row = 0;
    for (; row < root->rowCount(); row++, (*index)++)
    {
        QStandardItem* it = root->child(row);

        if (item->text().compare(it->text()) == 0)
        {
            return true;
        }

        if (it->hasChildren())
        {
            if (this->findIndex(item, it, index)) return true;
        }
    }
    (*index)--;
    return false;
}

void ECFParameterTreeWidget::onItemChanged(QStandardItem *item)
{
    int id = item->data(ECFParameterTreeDelegate::EditorRole).toInt();
    if (item->isCheckable()
            ||  this->m_values[id]->metadata.datatype.at(0) == ECFDataType::OPTION)
    {
        emit itemChanged();
    }

    this->resizeCellsToContent();
}

bool ECFParameterTreeWidget::isValid()
{
    bool valid = true;
    auto validateItems = [&] (QStandardItem* item) {
        bool ok;
        int id = item->data(ECFParameterTreeDelegate::EditorRole).toInt(&ok);
        if (ok && !item->isCheckable() && this->m_values[id]->metadata.ismandatory())
        {
            IValidatableObject *validatable;
            if ((validatable = dynamic_cast<IValidatableObject*>(this->m_delegate->editorFactory(id))))
            {
                valid = validatable->isValid();
                this->m_err = QString::fromStdString(this->m_values[id]->metadata.description[0])
                        + validatable->errorMessage();
            }
            else if (item->text().isEmpty())
            {
                valid = false;
                this->m_err = tr("Empty property %1")
                        .arg(QString::fromStdString(this->m_values[id]->metadata.description[0]));
            }
            if (!valid) return false;
        }
        return true;
    };
    this->applyForAllItems(this->m_root, validateItems);

    return valid;
}

QString ECFParameterTreeWidget::errorMessage()
{
    return this->m_err;
}

void ECFParameterTreeWidget::applyForAllItems(QStandardItem* root, std::function<bool(QStandardItem*)> apply)
{
    for (int row = 0; row < root->rowCount(); row++)
    {
        QStandardItem* item = root->child(row);

        if (item->hasChildren())
        {
            this->applyForAllItems(item, apply);
        }
        else {
            QStandardItem* value = root->child(row, 1);
            if (!apply(value)) return;
        }
    }
}

void ECFParameterTreeWidget::save()
{
    auto saveItems = [&] (QStandardItem* item) {
        bool ok;
        int id = item->data(ECFParameterTreeDelegate::EditorRole).toInt(&ok);
        if (item->isCheckable())
        {
            if (item->checkState() == Qt::Checked)
                this->m_values[id]->setValue("TRUE");
            else
                this->m_values[id]->setValue("FALSE");
        }
        else if (ok)
        {
            this->m_values[id]
                    ->setValue(item->text()
                               .toStdString());
        }
        return true;
    };
    this->applyForAllItems(this->m_root, saveItems);
}

void ECFParameterTreeWidget::saveState()
{
    this->m_stored_values.clear();

    for (int i = 0; i < this->m_values.size(); i++)
    {
        this->m_stored_values << this->m_values[i]->getValue();
    }
}

void ECFParameterTreeWidget::restoreState()
{
    if (!this->m_stored_values.size()) return;

    auto restore = [&] (QStandardItem* item)
    {
        bool ok;
        int i = item->data(ECFParameterTreeDelegate::EditorRole).toInt(&ok);
        if (!ok) return true;
        this->m_values[i]->setValue(this->m_stored_values[i]);
//        if (item->isCheckable())
//        {
//            if (this->m_values[i]->getValue().compare("TRUE") == 0)
//                item->setCheckState(Qt::Checked);
//            else
//                item->setCheckState(Qt::Unchecked);
//        }
//        else
//        {
//            item->setText(QString::fromStdString(this->m_values[i]->getValue()));
//        }
        return true;
    };

    this->applyForAllItems(m_root, restore);
}
