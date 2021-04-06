#ifndef WORKSPACEWINDOW_H
#define WORKSPACEWINDOW_H

#include <config/ecf/ecf.h>
#include <QMainWindow>

#include "config/configuration.h"
#include "input/sortedinput.h"
#include "declarations/datasetswidget.h"
#include "mesh/meshwidget.h"
#include "mesh/regionpickerwidget.h"
#include "workflowwidget.h"

#include "parallel/mpimanager.h"

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

    ECF* m_ecf;
    Mesh* m_mesh;

    WorkflowWidget* m_workflow;

    DataSetsWidget* m_datasets = nullptr;
    MeshWidget* m_mesh3D = nullptr;
    RegionPickerWidget* m_regions = nullptr;

    void initUi();
    void initPanels();
};

}

#endif // WORKSPACEWINDOW_H
