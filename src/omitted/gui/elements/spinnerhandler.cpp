#include "spinnerhandler.h"
#include "ui_spinnerhandler.h"

using namespace espreso;

SpinnerHandler::SpinnerHandler(ECFParameter* data,
                               bool withLabel,
                               QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SpinnerHandler)
{   
    this->m_data = data;

    ui->setupUi(this);

    ui->label->setText(QString::fromStdString(data->metadata.description[0]));
    if (!withLabel) ui->label->hide();

    ui->spinBox->setValue(QString::fromStdString(data->getValue()).toInt());
    connect(ui->spinBox, SIGNAL(valueChanged(int)),
            this, SLOT(onSpinnerValueChanged(int)));
}

SpinnerHandler::~SpinnerHandler()
{
    delete ui;
}

void SpinnerHandler::setValue(const QString& val)
{
    ui->spinBox->setValue(val.toInt());
}

void SpinnerHandler::save()
{
    this->m_data->setValue(QString::number(ui->spinBox->value()).toStdString());
}

void SpinnerHandler::saveState()
{
    this->m_saved_state = this->m_data->getValue();
}
void SpinnerHandler::restoreState()
{
    this->m_data->setValue(this->m_saved_state);
    ui->spinBox->setValue(QString::fromStdString(this->m_saved_state).toInt());
}

bool SpinnerHandler::isValid()
{
    return true;
}

QString SpinnerHandler::errorMessage()
{
    return tr("%1: Invalid value.").arg(ui->label->text());
}

void SpinnerHandler::onSpinnerValueChanged(int val)
{
    emit valueChanged(val);
}
