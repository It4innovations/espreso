
#ifndef LOADSTEPWIDGET_H
#define LOADSTEPWIDGET_H

#include "gui/elements/scrollecfobjectwidget.h"
#include "gui/elements/fieldhandler.h"
#include "gui/elements/regionpropertywidget.h"

namespace espreso
{

class LoadstepWidget : public ScrollECFObjectWidget
{
    Q_OBJECT

public:
    explicit LoadstepWidget(size_t id, QWidget* parent = 0);

protected:
    void drawObject(ECFObject*, int = 0) override;
    virtual void performBeforeRedraw() override;

private:
    ECFParameter* m_loadstep;

    RegionPropertyWidget* m_properties;
};

}

#endif // LOADSTEPWIDGET_H
