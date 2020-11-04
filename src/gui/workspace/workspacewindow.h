
#ifndef WORKSPACEWINDOW_H
#define WORKSPACEWINDOW_H

#include "workflowwidget.h"
#include "gui/declarations/datasetswidget.h"
#include "gui/mesh/meshwidget.h"
#include "gui/mesh/regionpickerwidget.h"
#include "gui/parallel/mpimanager.h"

#include "config/configuration.h"

#include <QMainWindow>

namespace espreso
{

namespace Ui {
class WorkspaceWindow;
}

class WorkspaceWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit WorkspaceWindow(MpiManager* manager, QWidget *parent = 0);
    ~WorkspaceWindow();

    void init();

private slots:
    void onPhysicsChanged(ECFObject* physics);
    void onInputChanged();
    void on_btnOpen_pressed();
    void on_btnSave_pressed();

private:
    Ui::WorkspaceWindow *ui;

    MpiManager* m_manager;

    WorkflowWidget* m_workflow;

    DataSetsWidget* m_datasets = nullptr;
    MeshWidget* m_mesh3D = nullptr;
    RegionPickerWidget* m_regions = nullptr;

    void initUi();
    void initPanels();
};

}

#endif // WORKSPACEWINDOW_H
