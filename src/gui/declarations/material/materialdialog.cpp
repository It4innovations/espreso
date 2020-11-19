
#include "materialdialog.h"
#include "materialpropertytablewidget.h"
#include "ui_materialdialog.h"

#include "gui/elements/optionhandler.h"
#include "gui/elements/formwidget.h"

#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>

using namespace espreso;

MaterialDialog::MaterialDialog(ECFObject* material,
                               const QVector<std::string>& materialNames,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MaterialDialog)
{
    ui->setupUi(this);

    this->m_widget = new MaterialWidget(material, materialNames, this);
    ui->frameLayout->addWidget(m_widget);
    this->m_widget->init();
}

MaterialDialog::~MaterialDialog()
{
    delete ui;
}

void MaterialDialog::accept()
{
    if (!this->m_widget->isValid()) return;

    QDialog::accept();
}
