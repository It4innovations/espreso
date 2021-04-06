#ifndef MATERIALWIDGET_H
#define MATERIALWIDGET_H

#include "elements/scrollecfobjectwidget.h"
#include "elements/formwidget.h"

namespace espreso
{

class MaterialWidget : public ScrollECFObjectWidget
{
    Q_OBJECT
public:
    MaterialWidget(MaterialConfiguration* material,
                   const QVector<std::string>& materialNames,
                   QWidget *parent = 0);

    bool isValid() override;


private:
    QVector<std::string> m_names;
    QWidget* m_widget;
};

}
#endif // MATERIALWIDGET_H
