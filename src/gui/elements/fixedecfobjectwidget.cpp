
#include "fixedecfobjectwidget.h"

using namespace espreso;

FixedECFObjectWidget::FixedECFObjectWidget(ECFObject* obj, QWidget* parent) :
    ECFObjectWidget(obj, parent)
{
}

QWidget* FixedECFObjectWidget::initContainer()
{
    return this->m_widget;
}
