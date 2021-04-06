#ifndef LOADSTEPWIDGET_H
#define LOADSTEPWIDGET_H

#include "elements/scrollecfobjectwidget.h"
#include "elements/fieldhandler.h"
#include "elements/regionpropertywidget.h"

namespace espreso
{

class LoadstepWidget : public ScrollECFObjectWidget
{
    Q_OBJECT

public:
    explicit LoadstepWidget(size_t id, Mesh* mesh, ECFObject* physics, QWidget* parent = 0);

protected:
    void drawObject(ECFObject*, int = 0) override;
    virtual void performBeforeRedraw() override;

private:
    ECFObject* m_physics;
    ECFParameter* m_loadstep;
    Mesh* m_mesh;

    RegionPropertyWidget* m_properties;
};

}

#endif // LOADSTEPWIDGET_H
