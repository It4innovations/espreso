
#include "loadstepwidget.h"

#include "gui/validators/validatorfactory.h"
#include "gui/elements/spinnerhandler.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

using namespace espreso;

LoadstepWidget::LoadstepWidget(size_t id, QWidget* parent) :
    ScrollECFObjectWidget(info::ecf->getPhysics()->ecfdescription, parent)
{
    auto *ph = info::ecf->getPhysics();
    this->m_loadstep = ph->ecfdescription->getParameter("load_steps_settings")->getParameter(QString::number(id).toStdString());
    this->m_obj = static_cast<ECFObject*>(m_loadstep);
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
            this->m_properties = new RegionPropertyWidget(tr("Boundary conditions"),
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
