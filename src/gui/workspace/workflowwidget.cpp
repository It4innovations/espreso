
#include "workflowwidget.h"
#include "ui_workflowwidget.h"

#include "loadstepwidget.h"
#include "regionmaterialswidget.h"
#include "outputconfigurationwidget.h"

#include "esinfo/ecfinfo.h"

#include <QFileDialog>
#include <QLabel>
#include <QComboBox>
#include <QDebug>
#include <QScrollArea>
#include <QMessageBox>

using namespace espreso;

WorkflowWidget::WorkflowWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WorkflowWidget)
{
    ui->setupUi(this);
}

WorkflowWidget::~WorkflowWidget()
{
    delete ui;
}

void WorkflowWidget::setData()
{   
    int tabs = ui->workflow->count();
    for (int i = 1; i < tabs; i++)
    {
        ui->workflow->removeTab(1);
    }

    this->m_physicsTab = nullptr;

    this->createInput();

    this->createPhysicsTab();

    this->createMaterialsTab();

    this->createLoadstepsTabs();

    this->m_loadsteps_fst_tab_index = ui->workflow->count();

    this->createOutputTab();
}

void WorkflowWidget::createInput()
{
    ECFParameter* i = info::ecf->ecfdescription->getParameter(&info::ecf->input_type);
    if (!this->m_inputBox_filled)
    {
        for (auto opt = i->metadata.options.begin(); opt != i->metadata.options.end(); ++opt)
        {
            ui->cmbInput->addItem( QString::fromStdString((*opt).description) );
        }
        connect(ui->cmbInput, SIGNAL(currentIndexChanged(int)),
                this, SLOT(onInputChange(int)));
        this->m_inputBox_filled = true;
    }
    int index = 0;
    for (auto opt = i->metadata.options.begin(); opt != i->metadata.options.end(); ++opt)
    {
        if ((*opt).name.compare(i->getValue()) == 0) break;
        index++;
    }
    ui->cmbInput->setCurrentIndex(index);

    this->onInputChange(ui->cmbInput->currentIndex());
}

void WorkflowWidget::onInputChange(int index)
{
    if (this->m_inputWidget != nullptr)
    {
        this->m_inputWidget->hide();
        ui->inputConfigLayout->removeWidget(this->m_inputWidget);
    }

    ECFParameter* input = info::ecf->ecfdescription->getParameter(&info::ecf->input_type);
    input->setValue(input->metadata.options[index].name);

    this->m_inputWidget = new InputWidget(this->input());
    this->m_inputWidget->init();
    ui->inputConfigLayout->addWidget(this->m_inputWidget);
}

ECFObject* WorkflowWidget::input()
{
    return info::ecf->input.ecfdescription;
}

void WorkflowWidget::createPhysicsTab()
{
    PhysicsWidget* pw = new PhysicsWidget(this);
    pw->init();

    this->m_phyDetail = pw;
    this->m_loadsteps = QString::fromStdString(
                pw->activePhysics()
                        ->getParameter("load_steps")
                        ->getValue()
                ).toInt();
    connect(pw, SIGNAL(loadstepsChanged(int)), this, SLOT(onLoadstepsChange(int)));
    connect(pw, SIGNAL(physicsChanged(ECFObject*)), this, SLOT(onPhysicsChange(ECFObject*)));

    ui->workflow->addTab(pw, QLatin1String("Physics"));
    ui->workflow->setCurrentIndex(ui->workflow->count() - 1);
}

void WorkflowWidget::createMaterialsTab()
{
    RegionMaterialsWidget* rmw = new RegionMaterialsWidget(this);
    ui->workflow->addTab(rmw, QLatin1String("Materials"));
}

void WorkflowWidget::createOutputTab()
{
    OutputConfigurationWidget* ocw = new OutputConfigurationWidget(info::ecf->output.ecfdescription, this);
    ocw->init();
    ui->workflow->addTab(ocw, QLatin1String("Output"));
}

void WorkflowWidget::createLoadstepsTabs()
{
    for (int i = 0; i < m_loadsteps; i++)
    {
        LoadstepWidget* lsw = new LoadstepWidget(i + 1, this);
        lsw->init();
        ui->workflow->addTab(lsw, tr("Loadstep %1").arg(i + 1));
    }
}

void WorkflowWidget::onLoadstepsChange(int loadsteps)
{
    int delta = this->m_loadsteps - loadsteps;
    if (delta > 0)
    {
        //DELETE LOADSTEPS

        ui->workflow->removeTab(this->m_loadsteps + this->m_loadsteps_fst_tab_index - 2);

        this->m_loadsteps--;
    }
    else if (delta < 0)
    {
        //ADD LOADSTEPS

        LoadstepWidget* lsw = new LoadstepWidget(++this->m_loadsteps, this);
        lsw->init();
        ui->workflow->insertTab(this->m_loadsteps + this->m_loadsteps_fst_tab_index - 2,
                                lsw,
                                tr("Loadstep %1").arg(this->m_loadsteps));
    }
}

PhysicsConfiguration* WorkflowWidget::activePhysics()
{
    return info::ecf->getPhysics();
}

void WorkflowWidget::onPhysicsChange(ECFObject *physics)
{

    int tabs = ui->workflow->count();
    for (int i = 2; i < tabs; i++)
    {
        ui->workflow->removeTab(2);
    }

    this->createMaterialsTab();

    this->m_loadsteps = QString::fromStdString(
                this->m_phyDetail
                        ->activePhysics()
                            ->getParameter("load_steps")
                                ->getValue()
                ).toInt();

    this->createLoadstepsTabs();

    this->m_loadsteps_fst_tab_index = ui->workflow->count();

    this->createOutputTab();

    emit physicsChanged(physics);
}

void WorkflowWidget::on_btnLoad_pressed()
{
    if (!this->m_inputWidget->isValid())
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(this->m_inputWidget->errorMessage());
        msg.exec();
        return;
    }

    this->m_inputWidget->save();
    emit inputChanged();
}

void WorkflowWidget::save()
{
   this->m_inputWidget->save();

   for (int i = 1; i < ui->workflow->count(); i++)
   {
       QWidget* widget = ui->workflow->widget(i);
       RegionMaterialsWidget* rmw = dynamic_cast<RegionMaterialsWidget*>(widget);
       if (rmw == nullptr)
       {
            static_cast<ECFObjectWidget*>(widget)->save();
       }
   }
}

bool WorkflowWidget::isValid()
{
    if (!this->m_inputWidget->isValid())
    {
        this->m_errmsg = this->m_inputWidget->errorMessage();
        return false;
    }

    for (int i = 1; i < ui->workflow->count(); i++)
    {
        QWidget* widget = ui->workflow->widget(i);
        RegionMaterialsWidget* rmw = dynamic_cast<RegionMaterialsWidget*>(widget);
        if (rmw == nullptr)
        {
             ECFObjectWidget *ecfwidget = static_cast<ECFObjectWidget*>(widget);
             if (!ecfwidget->isValid())
             {
                 this->m_errmsg = ecfwidget->errorMessage();
                 return false;
             }
        }
    }

    this->m_errmsg = "";
    return true;
}

QString WorkflowWidget::errorMessage()
{
    return this->m_errmsg;
}
