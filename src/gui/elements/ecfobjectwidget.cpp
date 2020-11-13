
#include "ecfobjectwidget.h"
#include "ui_ecfobjectwidget.h"

#include "gui/declarations/datatypeeditwidget.h"

#include <QMessageBox>
#include <QScrollArea>
#include <QLabel>

using namespace espreso;

ECFObjectWidget::ECFObjectWidget(ECFObject* object, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ECFObjectWidget)
{
    this->m_obj = object;
    ui->setupUi(this);
    this->m_parameter_groups.push(0);
}

ECFObjectWidget::~ECFObjectWidget()
{
    delete ui;
}

void ECFObjectWidget::init()
{
    this->drawMe();
}

void ECFObjectWidget::drawMe()
{
    this->m_savables.clear();
    this->m_validatables.clear();
    this->m_parameter_tree = nullptr;
    this->m_parameter_groups.clear();

    this->m_widget = this->initContainerWidget();
    this->m_container = this->initContainer();
    ui->layout->addWidget(this->m_container);

    this->drawObject(this->m_obj);
}

QWidget* ECFObjectWidget::initContainerWidget()
{
    QWidget* widget = new QWidget;
    QVBoxLayout* w_layout = new QVBoxLayout;
    w_layout->setMargin(0);
    w_layout->setSpacing(0);
    widget->setLayout(w_layout);

    return widget;
}

void ECFObjectWidget::drawObject(ECFObject *obj, int parentGroupId)
{
    this->m_parameter_tree = this->createParameterWidget(this->m_widget, this->m_parameter_tree);
    this->m_parameter_groups.push(
                this->m_parameter_tree->createGroup(
                    QString::fromStdString(obj->metadata.description.at(0)),
                    parentGroupId
                    ));

    this->processParameters(obj, this->m_widget);
}

OptionHandler* ECFObjectWidget::createOption(ECFParameter* option, QWidget* parent, bool withLabel)
{
    OptionHandler* handler = new OptionHandler(option, parent, withLabel);
    connect(handler, SIGNAL(optionChanged()), this, SLOT(redraw()));

    return handler;
}

BoolHandler* ECFObjectWidget::createBool(ECFParameter* param, QWidget* parent)
{
    BoolHandler* handler = new BoolHandler(param, parent);
    connect(handler, SIGNAL(stateChanged()), this, SLOT(redraw()));

    return handler;
}

void ECFObjectWidget::createHeadline(ECFObject* obj, QWidget* widget)
{
    if (!obj->metadata.description.size())
        return;

    QString lblNameText = QString::fromStdString(obj->metadata.description.at(0));
    QLabel* lblName = new QLabel(lblNameText);
    widget->layout()->addWidget(lblName);
}

void ECFObjectWidget::redraw()
{
    this->performBeforeRedraw();

    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    this->m_container->hide();
    ui->layout->removeWidget(this->m_container);

    this->drawMe();
}

bool ECFObjectWidget::validate()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            QMessageBox msg;
            msg.setWindowTitle(tr("Error"));
            msg.setText(obj->errorMessage());
            msg.exec();
            return false;
        }
    }

    return true;
}

bool ECFObjectWidget::isValid()
{
    foreach (IValidatableObject* obj, m_validatables) {
        if (!obj->isValid())
        {
            this->m_errormsg = obj->errorMessage();
            return false;
        }
    }

    return true;
}

QString ECFObjectWidget::errorMessage()
{
    return this->m_errormsg;
}

void ECFObjectWidget::save()
{
    foreach (ISavableObject* obj, m_savables)
    {
        obj->save();
    }
}

void ECFObjectWidget::setDrawHeadline(bool draw)
{
    this->m_draw_headlines = draw;
}

void ECFObjectWidget::processParameters(ECFObject* obj, QWidget* widget)
{
    int parentGroupId = this->m_parameter_groups.pop();

    // Iterating object parameters and creating UI elements for them
    for (auto parameter = obj->parameters.cbegin();
         parameter != obj->parameters.cend();
         ++parameter)
    {
        if (!(*parameter)->metadata.isallowed())
            continue;

        if ( (*parameter)->isObject() )
        {
            this->drawObject(static_cast<ECFObject*>(*parameter), parentGroupId);
        }
        else if ((*parameter)->metadata.datatype.size() == 1)
        {
            ECFDataType type = (*parameter)->metadata.datatype.at(0);

            if ( type == ECFDataType::OPTION
                 || type == ECFDataType::ENUM_FLAGS )
            {
                this->m_parameter_tree = this->processOptionEnum(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::BOOL)
            {
                this->m_parameter_tree = this->processBool(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::EXPRESSION)
            {
                this->m_parameter_tree = this->processExpression(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::STRING)
            {
                this->m_parameter_tree = this->processString(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::FLOAT)
            {
                this->m_parameter_tree = this->processFloat(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::NONNEGATIVE_INTEGER)
            {
                this->m_parameter_tree = this->processNonnegativeInteger(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::POSITIVE_INTEGER)
            {
                this->m_parameter_tree = this->processPositiveInteger(*parameter, m_parameter_tree, widget, parentGroupId);
            }
            else if (type == ECFDataType::BOUNDARY_REGION)
            {
                this->m_parameter_tree = this->processRegion(*parameter, m_parameter_tree, widget, parentGroupId);
            }
        }
    }

    if (m_parameter_tree != nullptr) m_parameter_tree->resizeCellsToContent();
}


ECFParameterTreeWidget* ECFObjectWidget::createParameterWidget(QWidget* container, ECFParameterTreeWidget* tree)
{
    if (tree != nullptr) return tree;

    ECFParameterTreeWidget* treeWidget = new ECFParameterTreeWidget;
    this->m_savables.append(treeWidget);
    this->m_validatables.append(treeWidget);
    container->layout()->addWidget(treeWidget);
    connect(treeWidget, SIGNAL(itemChanged()), this, SLOT(redraw()));

    return treeWidget;
}

ECFParameterTreeWidget* ECFObjectWidget::processParameter(ECFParameter *param, ECFParameterTreeWidget *tree, QWidget *parent, int groupId)
{
    ECFParameterTreeWidget* treeWidget = this->createParameterWidget(parent, tree);
    treeWidget->add(static_cast<ECFValue*>(param), groupId);

    return treeWidget;
}

ECFParameterTreeWidget* ECFObjectWidget::processExpression(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    ECFParameterTreeWidget* treeWidget = this->createParameterWidget(widget, tree);
    DataTypeEditWidgetFactory *factory = new DataTypeEditWidgetFactory(parameter);
    treeWidget->addWithEditorFactory(static_cast<ECFValue*>(parameter), factory, groupId);

    return treeWidget;
}

ECFParameterTreeWidget* ECFObjectWidget::processOptionEnum(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processBool(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processString(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processFloat(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processNonnegativeInteger(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processPositiveInteger(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}

ECFParameterTreeWidget* ECFObjectWidget::processRegion(ECFParameter* parameter, ECFParameterTreeWidget* tree, QWidget* widget, int groupId)
{
    return this->processParameter(parameter, tree, widget, groupId);
}
