
#ifndef ECFOBJECTWIDGET_H
#define ECFOBJECTWIDGET_H

#include "isavableobject.h"
#include "ivalidatableobject.h"
#include "optionhandler.h"
#include "boolhandler.h"
#include "formwidget.h"
#include "regionpropertywidget.h"
#include "ecfvaluetablewidget.h"
#include "ecfparametertreewidget.h"
#include "gui/declarations/material/materialpropertytablewidget.h"

#include "config/ecf/ecf.h"

#include <QWidget>
#include <QStack>

namespace espreso
{

namespace Ui {
class ECFObjectWidget;
}

class ECFObjectWidget : public QWidget,
        public IValidatableObject, public ISavableObject
{
    Q_OBJECT

public:
    explicit ECFObjectWidget(ECFObject* object, QWidget *parent = 0);
    virtual ~ECFObjectWidget();
    void init();

    virtual bool isValid() override;
    virtual QString errorMessage() override;

    virtual void save() override;

    void setDrawHeadline(bool draw);

protected slots:
    void redraw();

protected:
    Ui::ECFObjectWidget *ui;

    ECFObject* m_obj;

    QWidget* m_container;
    QWidget* m_widget;

    QVector<ISavableObject*> m_savables;
    QVector<IValidatableObject*> m_validatables;
    bool validate();

    bool m_draw_headlines = true;

    QString m_errormsg;

    virtual QWidget* initContainerWidget();
    virtual QWidget* initContainer() = 0;
    virtual void performBeforeRedraw() = 0;

    virtual void drawObject(ECFObject* obj, int parentGroupId = 0);
    void drawMe();

    ECFParameterTreeWidget* m_parameter_tree = nullptr;
    QStack<int> m_parameter_groups;

    void processParameters(ECFObject*, QWidget*);

    virtual ECFParameterTreeWidget* processParameter(ECFParameter* param,
                                                     ECFParameterTreeWidget* table,
                                                     QWidget* parent,
                                                     int groupId = 0);

    virtual ECFParameterTreeWidget* processExpression(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processOptionEnum(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processBool(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processString(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processFloat(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processNonnegativeInteger(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processPositiveInteger(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);
    virtual ECFParameterTreeWidget* processRegion(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0);

    ECFParameterTreeWidget* createParameterWidget(QWidget*, ECFParameterTreeWidget*);

    OptionHandler* createOption(ECFParameter*, QWidget* = 0, bool = true);
    BoolHandler* createBool(ECFParameter*, QWidget* = 0);
    void createHeadline(ECFObject*, QWidget*);
};

}

#endif // ECFOBJECTWIDGET_H
