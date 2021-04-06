#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include "elements/fixedecfobjectwidget.h"

namespace espreso
{

class InputWidget : public FixedECFObjectWidget
{
    Q_OBJECT
public:
    InputWidget(ECFObject* obj, QWidget* parent = 0);

protected:
    void drawObject(ECFObject*, int = 0) override;
//    virtual ECFValueTableWidget* processString(ECFParameter*, ECFValueTableWidget*, QWidget*) override;
    virtual ECFParameterTreeWidget* processString(ECFParameter*, ECFParameterTreeWidget*, QWidget*, int = 0) override;
};

}

#endif // INPUTWIDGET_H
