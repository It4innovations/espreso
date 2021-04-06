#include "physicswidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "validators/validatorfactory.h"
#include "elements/spinnerhandler.h"

using namespace espreso;

PhysicsWidget::PhysicsWidget(ECF* ecf, Mesh* mesh, QWidget* parent) :
    ScrollECFObjectWidget(ecf, parent)
{
    this->m_ecf = ecf;
    this->m_mesh = mesh;
}

QWidget* PhysicsWidget::initContainer()
{
    ECFParameter* physics = m_ecf->getParameter("physics");

    QWidget* cmbLine = new QWidget(this->m_widget);
    QHBoxLayout* layout = new QHBoxLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    cmbLine->setLayout(layout);
    this->m_widget->layout()->addWidget(cmbLine);

    QComboBox* cmbPhysics = new QComboBox(cmbLine);
    QLabel* lblPhysics = new QLabel(tr("Physics:"), cmbLine);
    lblPhysics->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
    layout->addWidget(lblPhysics);
    layout->addWidget(cmbPhysics);

    int active = 0;
    int index = 0;
    for (auto option = physics->metadata.options.begin();
         option != physics->metadata.options.end();
         ++option, ++index)
    {
        QString name = QString::fromStdString(option->name);
        cmbPhysics->addItem(name);

        if (option->name.compare(physics->getValue()) == 0) active = index;
    }

    cmbPhysics->setCurrentIndex(active);
    this->m_physics = cmbPhysics;
    this->m_obj = this->physics(active);
    connect(cmbPhysics, SIGNAL(currentIndexChanged(int)),
            this, SLOT(onPhysicsChange(int)));

    return ScrollECFObjectWidget::initContainer();
}

ECFObject* PhysicsWidget::physics(int index)
{
    switch (index)
    {
        case 0:
            return &m_ecf->heat_transfer_2d;
        case 1:
            return &m_ecf->heat_transfer_3d;
        case 2:
            return &m_ecf->structural_mechanics_2d;
        case 3:
            return &m_ecf->structural_mechanics_3d;
        default:
            qCritical() << tr("WorkflowWidget: Invalid index of physics!");
            return nullptr;
    }
}

void PhysicsWidget::onPhysicsChange(int index)
{
    ECFParameter* physics = m_ecf->getParameter("physics");
    physics->setValue(physics->metadata.options[index].name);

    this->m_obj = this->physics(index);
    this->m_properties = nullptr;

    this->redraw();

    emit physicsChanged(this->m_obj);
}

void PhysicsWidget::drawObject(ECFObject* obj, int parentGroupId)
{
    if (obj->name.compare("material_set") == 0
            || obj->name.compare("load_steps_settings") == 0
            || obj->name.compare("materials") == 0)
        return;

    if ( obj->metadata.datatype.size() == 2 )
    {
        if (this->m_properties == nullptr)
        {
            this->m_properties =
                    new RegionPropertyWidget(m_mesh,
                                             static_cast<PhysicsConfiguration*>(this->activePhysics()),
                                             tr("Element region properties"),
                                             this->m_container);
        }
        this->m_properties->addProperty(obj);
        this->m_widget->layout()->addWidget(m_properties);

        return;
    }

    ScrollECFObjectWidget::drawObject(obj, parentGroupId);
}

void PhysicsWidget::performBeforeRedraw()
{
    this->m_properties = nullptr;
}

ECFParameterTreeWidget* PhysicsWidget::processPositiveInteger(ECFParameter* parameter,
                                                              ECFParameterTreeWidget* table,
                                                              QWidget* widget,
                                                              int groupId)
{
    ECFParameterTreeWidget* tw = this->createParameterWidget(widget, table);
    if (parameter->name.compare("load_steps") == 0)
    {
        SpinnerHandler* handler = new SpinnerHandler(parameter, false, widget);
        connect(handler, SIGNAL(valueChanged(int)), this, SLOT(onLoadstepsChange(int)));
        tw->addWithWidget(static_cast<ECFValue*>(parameter), handler, groupId);
        this->m_savables.append(handler);
        this->m_validatables.append(handler);

        return tw;
    }
    else
    {
        return ECFObjectWidget::processPositiveInteger(parameter, table, widget, groupId);
    }
}

void PhysicsWidget::onLoadstepsChange(int loadsteps)
{
    emit this->loadstepsChanged(loadsteps);
}

ECFObject* PhysicsWidget::activePhysics()
{
    return this->physics(m_physics->currentIndex());
}
