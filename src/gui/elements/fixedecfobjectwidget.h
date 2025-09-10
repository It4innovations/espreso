
#ifndef FIXEDECFOBJECTWIDGET_H
#define FIXEDECFOBJECTWIDGET_H

#include "ecfobjectwidget.h"

namespace espreso
{

class FixedECFObjectWidget : public ECFObjectWidget
{
    Q_OBJECT

public:
    FixedECFObjectWidget(ECFObject* obj, QWidget* parent = 0);

protected:
    QWidget* initContainer() override;
    void performBeforeRedraw() override {}
};

}

#endif // FIXEDECFOBJECTWIDGET_H
