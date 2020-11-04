
#include "fieldhandler.h"
#include "ui_fieldhandler.h"

using namespace espreso;

FieldHandler::FieldHandler(ECFParameter* data,
                           const ValidatorFactory* validator,
                           bool withLabel,
                           QWidget *parent) :
    QWidget(parent),
    ui(new Ui::FieldHandler)
{
    this->m_data = data;

    ui->setupUi(this);

    ui->label->setText(QString::fromStdString(data->metadata.description[0]));
    if (!withLabel) ui->label->hide();

    if (validator != nullptr)
        ui->field->setValidator(validator->create(this));

    ui->field->setText(QString::fromStdString(data->getValue()));
}

FieldHandler::~FieldHandler()
{
    delete ui;
}

void FieldHandler::setValue(const QString& val)
{
    ui->field->setText(val);
}

QString FieldHandler::value() const
{
    return ui->field->text();
}

void FieldHandler::save()
{
    this->m_data->setValue(ui->field->text().toStdString());
}

void FieldHandler::saveState()
{
    this->m_saved_state = this->m_data->getValue();
}
void FieldHandler::restoreState()
{
    this->m_data->setValue(this->m_saved_state);
    ui->field->setText(QString::fromStdString(this->m_saved_state));
}

bool FieldHandler::isValid()
{
    return true;
}

QString FieldHandler::errorMessage()
{
    return tr("%1: Invalid value.").arg(QString::fromStdString(m_data->name));
}
