#include "regionmaterialswidget.h"
#include "ui_regionmaterialswidget.h"

#include "elements/regionpropertywidget.h"

using namespace espreso;

RegionMaterialsWidget::RegionMaterialsWidget(Mesh* mesh, PhysicsConfiguration* physics, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::RegionMaterialsWidget)
{
    ui->setupUi(this);

    RegionPropertyWidget* rpw = new RegionPropertyWidget(mesh, physics, "", this);
    rpw->addProperty(static_cast<ECFObject*>(physics->getParameter("material_set")));

    ui->layout->addWidget(rpw);
}

RegionMaterialsWidget::~RegionMaterialsWidget()
{
    delete ui;
}
