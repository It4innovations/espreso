
#include "optionwidget.h"
#include "ui_optionwidget.h"

using namespace espreso;

OptionWidget::OptionWidget(const QStringList& options, QWidget *parent) :
    TextItemWidget(parent),
    ui(new Ui::OptionWidget)
{
    ui->setupUi(this);

    ui->cmb->addItems(options);

    this->m_options = options;

    connect(ui->cmb, SIGNAL(currentIndexChanged(QString)),
            this, SLOT(onIndexChanged(QString)));
}

OptionWidget::~OptionWidget()
{
    delete ui;
}

void OptionWidget::setText(const QString& text)
{
    int selected = this->m_options.indexOf(text);
    if (this->m_lastOption.isEmpty()) this->m_lastOption = text;
    ui->cmb->setCurrentIndex(selected);
}

QString OptionWidget::text()
{
    return ui->cmb->currentText();
}

void OptionWidget::onIndexChanged(QString newOption)
{
    if (this->m_lastOption.compare(newOption) != 0)
    {
        emit finished(this);
    }
    this->m_lastOption = newOption;
}
