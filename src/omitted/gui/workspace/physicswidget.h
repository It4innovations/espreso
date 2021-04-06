#ifndef PHYSICSWIDGET_H
#define PHYSICSWIDGET_H

#include <QPushButton>
#include <QComboBox>

#include "elements/scrollecfobjectwidget.h"
#include "elements/fieldhandler.h"
#include "elements/regionpropertywidget.h"

namespace espreso
{

class PhysicsWidget : public ScrollECFObjectWidget
{
    Q_OBJECT

public:
    explicit PhysicsWidget(ECF* ecf, Mesh* mesh, QWidget* parent = 0);

    ECFObject* activePhysics();

signals:
    void loadstepsChanged(int loadsteps);
    void physicsChanged(ECFObject* physics);

protected:
    QWidget* initContainer() override;
    virtual void drawObject(ECFObject* obj, int parentGroupId = 0) override;
    virtual void performBeforeRedraw() override;
    //ECFValueTableWidget* processPositiveInteger(ECFParameter *, ECFParameterTreeWidget *, QWidget *) override;
    ECFParameterTreeWidget* processPositiveInteger(ECFParameter *, ECFParameterTreeWidget *, QWidget *, int = 0) override;

private slots:
    void onPhysicsChange(int index);
    void onLoadstepsChange(int loadsteps);

private:
    ECF* m_ecf;
    Mesh* m_mesh;

    QComboBox* m_physics;

    RegionPropertyWidget* m_properties = nullptr;

    ECFObject* physics(int index);
};

}

#endif // PHYSICSWIDGET_H
