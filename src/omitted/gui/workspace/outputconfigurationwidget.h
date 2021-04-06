#ifndef OUTPUTCONFIGURATIONWIDGET_H
#define OUTPUTCONFIGURATIONWIDGET_H

#include "elements/scrollecfobjectwidget.h"

namespace espreso
{

class OutputConfigurationWidget : public ScrollECFObjectWidget
{
    Q_OBJECT

public:
    OutputConfigurationWidget(ECFObject* output, QWidget* parent = 0);

protected:
    virtual void drawObject(ECFObject*, int = 0) override;
};

}

#endif // OUTPUTCONFIGURATIONWIDGET_H
