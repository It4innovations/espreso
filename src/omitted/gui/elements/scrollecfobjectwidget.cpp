#include "scrollecfobjectwidget.h"

#include <QScrollArea>

using namespace espreso;

ScrollECFObjectWidget::ScrollECFObjectWidget(ECFObject* obj, QWidget* parent) :
    ECFObjectWidget(obj, parent)
{

}

QWidget* ScrollECFObjectWidget::initContainer()
{
    QScrollArea* area = new QScrollArea;

    area->setWidgetResizable(true);
    area->setWidget(this->m_widget);

    return area;
}
