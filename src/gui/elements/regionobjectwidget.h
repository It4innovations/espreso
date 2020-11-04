
#ifndef REGIONOBJECTWIDGET_H
#define REGIONOBJECTWIDGET_H

#include "fixedecfobjectwidget.h"

namespace espreso
{

class RegionObjectWidget : public FixedECFObjectWidget
{
    Q_OBJECT
public:
    explicit RegionObjectWidget(ECFObject* obj, QWidget *parent = nullptr);
};

}

#endif // REGIONOBJECTWIDGET_H
