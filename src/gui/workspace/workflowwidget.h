
#ifndef WORKFLOWWIDGET_H
#define WORKFLOWWIDGET_H

#include "physicswidget.h"
#include "inputwidget.h"

#include "config/configuration.h"

#include <QWidget>

namespace espreso
{

namespace Ui {
class WorkflowWidget;
}

class WorkflowWidget : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

public:
    explicit WorkflowWidget(QWidget *parent = 0);
    ~WorkflowWidget();

    void setData();
    PhysicsConfiguration* activePhysics();
    ECFObject* input(int index);

    virtual void save() override;
    virtual bool isValid() override;
    virtual QString errorMessage() override;

signals:
    void inputChanged();
    void physicsChanged(ECFObject* physics);

private slots:
    void onLoadstepsChange(int loadsteps);
    void onPhysicsChange(ECFObject* physics);
    void onInputChange(int index);

    void on_btnLoad_pressed();

private:
    Ui::WorkflowWidget *ui;

    bool m_inputBox_filled = false;
    void createInput();
    InputWidget* m_inputWidget = nullptr;

    void createPhysicsTab();
    QWidget* m_physicsTab;
    PhysicsWidget* m_phyDetail = nullptr;

    void createMaterialsTab();

    int m_loadsteps;
    int m_loadsteps_fst_tab_index;
    void createLoadstepsTabs();

    void createOutputTab();

    QString m_errmsg;
};

}

#endif // WORKFLOWWIDGET_H
