#include "boolhandler.h"
#include "ui_boolhandler.h"

#include <QDebug>

using namespace espreso;

BoolHandler::BoolHandler(ECFParameter* param, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BoolHandler)
{
    ui->setupUi(this);

    this->m_param = param;
//    QString label = QString::fromStdString(param->metadata.description[0]);
    ui->chck->setText("");

    if (param->getValue().compare("FALSE") == 0) ui->chck->setChecked(false);
    else ui->chck->setChecked(true);

    connect(ui->chck, SIGNAL(stateChanged(int)), this, SLOT(onStateChanged(int)));
}

BoolHandler::~BoolHandler()
{
    delete ui;
}

void BoolHandler::onStateChanged(int state)
{
    if (state) this->m_param->setValue("TRUE");
    else this->m_param->setValue("FALSE");

    emit stateChanged();
}
