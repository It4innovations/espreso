
#include "pathhandler.h"

#include <QHBoxLayout>
#include <QPushButton>
#include <QFileDialog>

using namespace espreso;

PathHandler::PathHandler(ECFParameter* path, QWidget* parent) :
    QWidget(parent)
{
    this->m_path = path;

    QHBoxLayout* layout = new QHBoxLayout;
    this->setLayout(layout);

    QPushButton* btn = new QPushButton(tr("Open file"), this);
    layout->addWidget(btn);
    connect(btn, SIGNAL(pressed()), this, SLOT(onPressed()));

    this->m_label = new QLabel(this);
    if (!path->getValue().empty())
        this->m_label->setText(QString::fromStdString(path->getValue()));
    layout->addWidget(this->m_label);
}

void PathHandler::onPressed()
{
    QString filename = QFileDialog::getOpenFileName(this,
                                                    tr("Open File"), ".");
    this->m_label->setText(filename);
    this->m_path->setValue(filename.toStdString());
}

bool PathHandler::isValid()
{
    if (this->m_path->getValue().empty())
        return false;

    return true;
}

QString PathHandler::errorMessage()
{
    return tr("\'%1\' is empty")
            .arg(QString::fromStdString(this->m_path->metadata.description[0]));
}
