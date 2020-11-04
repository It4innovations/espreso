
#ifndef REGIONMATERIALSWIDGET_H
#define REGIONMATERIALSWIDGET_H

#include "gui/elements/isavableobject.h"
#include "gui/elements/ivalidatableobject.h"

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
    explicit RegionMaterialsWidget(QWidget *parent = 0);
    ~RegionMaterialsWidget();

private:
    Ui::RegionMaterialsWidget *ui;
};

}

#endif // REGIONMATERIALSWIDGET_H
