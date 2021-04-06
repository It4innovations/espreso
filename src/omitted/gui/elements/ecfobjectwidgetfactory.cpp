#include "ecfobjectwidgetfactory.h"

using namespace espreso;

FixedECFObjectWidgetFactory::FixedECFObjectWidgetFactory(bool drawHeadlines) :
    ECFObjectWidgetFactory()
{
    this->m_draw_headlines = drawHeadlines;
}

ECFObjectWidget* FixedECFObjectWidgetFactory::create(ECFObject *object, QWidget *parent)
{
    FixedECFObjectWidget* w = new FixedECFObjectWidget(object, parent);
    w->setDrawHeadline(this->m_draw_headlines);
    return w;
}
