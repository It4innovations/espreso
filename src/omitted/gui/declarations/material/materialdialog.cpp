#include "materialdialog.h"
#include "ui_materialdialog.h"

#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>
#include "elements/optionhandler.h"
#include "elements/formwidget.h"
#include "materialpropertytablewidget.h"

using namespace espreso;

MaterialDialog::MaterialDialog(MaterialConfiguration* material,
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
