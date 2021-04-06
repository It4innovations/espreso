#include "outputconfigurationwidget.h"

#include <QVBoxLayout>
#include <QFormLayout>
#include <QScrollArea>
#include <QLabel>
#include <QDebug>

#include "elements/ecfobjectwidgetfactory.h"
#include "elements/integertabwidget.h"
#include <memory>

using namespace espreso;

OutputConfigurationWidget::OutputConfigurationWidget(ECFObject* output, QWidget* parent) :
    ScrollECFObjectWidget(output, parent)
{

}

void OutputConfigurationWidget::drawObject(ECFObject* obj, int groupId)
{
    if (obj->name.compare("monitoring") == 0)
    {
        std::unique_ptr<ECFObjectWidgetFactory> factory(new FixedECFObjectWidgetFactory(false));
        IntegerTabWidget* w = new IntegerTabWidget(obj, std::move(factory));
        this->m_widget->layout()->addWidget(w);
        this->m_savables.append(w);
        this->m_validatables.append(w);
        return;
    }

    ScrollECFObjectWidget::drawObject(obj, groupId);
}
