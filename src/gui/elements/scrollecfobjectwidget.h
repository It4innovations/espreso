
#ifndef SCROLLECFOBJECTWIDGET_H
#define SCROLLECFOBJECTWIDGET_H

#include "ecfobjectwidget.h"

namespace espreso
{

class ScrollECFObjectWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    ScrollECFObjectWidget(ECFObject* obj, QWidget* parent = 0);

protected:
    QWidget* initContainer() override;
    void performBeforeRedraw() override {}
};

}

#endif // SCROLLECFOBJECTWIDGET_H
