
#ifndef ECFOBJECTWIDGETFACTORY_H
#define ECFOBJECTWIDGETFACTORY_H

#include "ecfobjectwidget.h"
#include "fixedecfobjectwidget.h"

namespace espreso
{

struct ECFObjectWidgetFactory
{
public:
    virtual ~ECFObjectWidgetFactory() {}

    virtual ECFObjectWidget* create(ECFObject* object, QWidget *parent = 0) = 0;
protected:
    ECFObjectWidgetFactory() {}
};


class FixedECFObjectWidgetFactory : public ECFObjectWidgetFactory
{
public:
    FixedECFObjectWidgetFactory(bool drawHeadlines = true);
    ECFObjectWidget* create(ECFObject* object, QWidget *parent = 0) override;

private:
    bool m_draw_headlines;
};

}

#endif // ECFOBJECTWIDGETFACTORY_H
