
#ifndef PHYSICSWIDGET_H
#define PHYSICSWIDGET_H

#include "gui/elements/scrollecfobjectwidget.h"
#include "gui/elements/fieldhandler.h"
#include "gui/elements/regionpropertywidget.h"

#include "config/configuration.h"

#include <QPushButton>
#include <QComboBox>

namespace espreso
{

class PhysicsWidget : public ScrollECFObjectWidget
{
    Q_OBJECT

public:
    explicit PhysicsWidget(QWidget* parent = 0);

    ECFObject* activePhysics();

signals:
    void loadstepsChanged(int loadsteps);
    void physicsChanged(ECFObject* physics);

protected:
    QWidget* initContainer() override;
    void drawObject(ECFObject* obj, int parentGroupId = 0) override;
    void performBeforeRedraw() override;
    //ECFValueTableWidget* processPositiveInteger(ECFParameter *, ECFParameterTreeWidget *, QWidget *) override;
    ECFParameterTreeWidget* processPositiveInteger(ECFParameter *, ECFParameterTreeWidget *, QWidget *, int = 0) override;

private slots:
    void onPhysicsChange(int index);
    void onLoadstepsChange(int loadsteps);

private:
    QComboBox* m_physics;

    RegionPropertyWidget* m_properties = nullptr;

    ECFObject* physics(int index);
};

}

#endif // PHYSICSWIDGET_H
