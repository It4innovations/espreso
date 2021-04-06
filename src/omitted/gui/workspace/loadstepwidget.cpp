#include "loadstepwidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "validators/validatorfactory.h"
#include "elements/spinnerhandler.h"

using namespace espreso;

LoadstepWidget::LoadstepWidget(size_t id, Mesh* mesh, ECFObject* physics, QWidget* parent) :
    ScrollECFObjectWidget(physics, parent)
{
    this->m_physics = physics;
    this->m_loadstep = m_physics->getParameter("load_steps_settings")->getParameter(QString::number(id).toStdString());
    this->m_obj = static_cast<ECFObject*>(m_loadstep);
    this->m_mesh = mesh;
    this->m_properties = nullptr;
}

void LoadstepWidget::drawObject(ECFObject* obj, int groupId)
{
    if (obj->name.compare("material_set") == 0
            || obj->name.compare("load_steps_settings") == 0)
        return;

    if ( obj->metadata.datatype.size() == 2 || obj->metadata.description.size() == 2)
    {
        if (this->m_properties == nullptr)
        {
            this->m_properties = new RegionPropertyWidget(m_mesh,
                                                          static_cast<PhysicsConfiguration*>(m_physics),
                                                          tr("Boundary conditions"),
                                                          this->m_container);
        }
        this->m_properties->addProperty(obj);
        this->m_widget->layout()->addWidget(m_properties);
        return;
    }

    ScrollECFObjectWidget::drawObject(obj, groupId);
}

void LoadstepWidget::performBeforeRedraw()
{
    this->m_properties = nullptr;
}
