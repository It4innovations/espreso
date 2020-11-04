
#include "filepathwidget.h"
#include "ui_filepathwidget.h"

#include <QFileDialog>

using namespace espreso;

FilepathWidget::FilepathWidget(QWidget *parent) :
    TextItemWidget(parent),
    ui(new Ui::FilepathWidget)
{
    ui->setupUi(this);

    connect(ui->button, SIGNAL(pressed()), this, SLOT(onPressed()));
}

FilepathWidget::~FilepathWidget()
{
    delete ui;
}

void FilepathWidget::onPressed()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open File"), ".");
    if (filename.isEmpty()) return;

    ui->label->setText(filename);

    emit finished(this);
}

void FilepathWidget::setText(const QString& text)
{
    ui->label->setText(text);
}

QString FilepathWidget::text()
{
    return ui->label->text();
}
