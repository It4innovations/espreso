
#include "regionmaterialswidget.h"
#include "ui_regionmaterialswidget.h"

#include "gui/elements/regionpropertywidget.h"

#include "esinfo/ecfinfo.h"

using namespace espreso;

RegionMaterialsWidget::RegionMaterialsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::RegionMaterialsWidget)
{
    ui->setupUi(this);

    RegionPropertyWidget* rpw = new RegionPropertyWidget("", this);
    rpw->addProperty(static_cast<ECFObject*>(info::ecf->getPhysics()->ecfdescription->getParameter("material_set")));

    ui->layout->addWidget(rpw);
}

RegionMaterialsWidget::~RegionMaterialsWidget()
{
    delete ui;
}
