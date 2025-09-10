
#ifndef INTEGERTABWIDGET_H
#define INTEGERTABWIDGET_H

#include "ecfobjectwidgetfactory.h"
#include "isavableobject.h"
#include "ivalidatableobject.h"

#include "config/configuration.h"

#include <QWidget>
#include <QTabWidget>

namespace espreso
{

class IntegerTabWidget : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

public:
    IntegerTabWidget(ECFObject* map, std::unique_ptr<ECFObjectWidgetFactory> factory, QWidget* parent = 0);

    void save() override;
    bool isValid() override;
    QString errorMessage() override;

private slots:
    void onTabClosed(int index);
    void onAddPressed();

private:
    QTabWidget* m_tabwidget;

    ECFObject* m_map;
    std::unique_ptr<ECFObjectWidgetFactory> m_factory;

    QString m_errmsg;

    int m_key;

    void addParam(ECFParameter*);
};

}

#endif // INTEGERTABWIDGET_H
