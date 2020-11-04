
#include "workspacewindow.h"
#include "ui_workspacewindow.h"

#include "config/reader/reader.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"

#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>

#include <iostream>
#include <fstream>

using namespace espreso;

WorkspaceWindow::WorkspaceWindow(MpiManager* manager, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::WorkspaceWindow)
{
    this->m_manager = manager;
}

WorkspaceWindow::~WorkspaceWindow()
{
    delete ui;
}

void WorkspaceWindow::init()
{
    ui->setupUi(this);
    this->initUi();

    this->initPanels();
    this->m_workflow->setData();
}

void WorkspaceWindow::initUi()
{
    this->m_workflow = new WorkflowWidget(this);
    ui->right->layout()->addWidget(m_workflow);
    connect(m_workflow, SIGNAL(inputChanged()), this, SLOT(onInputChanged()));
    connect(m_workflow, SIGNAL(physicsChanged(ECFObject*)), this, SLOT(onPhysicsChanged(ECFObject*)));
    if (info::mpi::rank == 0) MeshWidget::initOGL();
}

void WorkspaceWindow::onInputChanged()
{
    this->initPanels();
    this->m_workflow->setData();
}

void WorkspaceWindow::initPanels()
{
    if (this->m_datasets != nullptr)
    {
        this->m_datasets->setVisible(false);
        ui->left->layout()->removeWidget(this->m_datasets);
        delete this->m_datasets;
        this->m_datasets = nullptr;
    }

    if (this->m_mesh3D != nullptr)
    {
        this->m_mesh3D->setVisible(false);
        ui->center->layout()->removeWidget(this->m_mesh3D);
        delete this->m_mesh3D;
        this->m_mesh3D = nullptr;
    }

    if (this->m_regions != nullptr)
    {
        this->m_regions->setVisible(false);
        ui->left->layout()->removeWidget(this->m_regions);
        delete this->m_regions;
        this->m_regions = nullptr;
    }

    ECFObject* materials = dynamic_cast<ECFObject*>(info::ecf->getPhysics()->ecfdescription->getParameter(&info::ecf->getPhysics()->materials));
    this->m_datasets = new DataSetsWidget(materials, tr("Data sets"), this);
    ui->left->layout()->addWidget(m_datasets);

    this->m_mesh3D = new MeshWidget(m_manager, this);
    ui->center->layout()->addWidget(m_mesh3D);
    m_mesh3D->setVisible(true);

    this->m_regions = new RegionPickerWidget(m_mesh3D, this);
    ui->left->layout()->addWidget(m_regions);

    // Set size of central widgets (data sets, 3D mesh, workflow)
    QList<int> sizes;
    sizes << 200 << 800 << 300;
    ui->splitter->setSizes(sizes);
}

void WorkspaceWindow::onPhysicsChanged(ECFObject *physics)
{
    ECFObject* materials = static_cast<ECFObject*>(physics->getParameter("materials"));
    this->m_datasets->setMaterials(materials);
}

void espreso::WorkspaceWindow::on_btnOpen_pressed()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open ESPRESO Configuration File"), ".", tr("ECF (*.ecf)"));
    if (filename.isEmpty()) return;

    this->m_manager->masterOpenECF(filename);

    this->initPanels();
    this->m_workflow->setData();
}

void espreso::WorkspaceWindow::on_btnSave_pressed()
{
    if (!this->m_workflow->isValid())
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(this->m_workflow->errorMessage());
        msg.exec();

        return;
    }

    this->m_workflow->save();

    QString path = QFileDialog::getSaveFileName(this, tr("Save Configuration As"), tr("espreso.ecf"));
    std::ofstream file;
    file.open(path.toStdString());

    ECFReader::store(*info::ecf->ecfdescription, file);
}
