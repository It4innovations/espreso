
#include "textwidget.h"
#include "ui_textwidget.h"

using namespace espreso;

TextWidget::TextWidget(QWidget *parent) :
    TextItemWidget(parent),
    ui(new Ui::TextWidget)
{
    ui->setupUi(this);
}

TextWidget::~TextWidget()
{
    delete ui;
}

void TextWidget::setValidator(QValidator* validator)
{
    ui->editor->setValidator(validator);
}

void TextWidget::setText(const QString& text)
{
    ui->editor->setText(text);
}

QString TextWidget::text()
{
    return ui->editor->text();
}
