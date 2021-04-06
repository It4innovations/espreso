#ifndef REGIONMATERIALSWIDGET_H
#define REGIONMATERIALSWIDGET_H

#include <config/ecf/ecf.h>
#include "config/configuration.h"
#include "mesh/mesh.h"

#include "elements/isavableobject.h"
#include "elements/ivalidatableobject.h"

#include <QWidget>

namespace espreso
{

namespace Ui {
class RegionMaterialsWidget;
}

class RegionMaterialsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit RegionMaterialsWidget(Mesh* mesh, PhysicsConfiguration* physics, QWidget *parent = 0);
    ~RegionMaterialsWidget();

private:
    Ui::RegionMaterialsWidget *ui;
};

}

#endif // REGIONMATERIALSWIDGET_H
