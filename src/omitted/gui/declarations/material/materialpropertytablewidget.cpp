#include "materialpropertytablewidget.h"
#include "ui_materialpropertytablewidget.h"

#include <QComboBox>
#include <QLabel>
#include <QLineEdit>

#include "validators/validatorfactory.h"
#include "tablewidget.h"

using namespace espreso;

MaterialPropertyTableWidget::MaterialPropertyTableWidget(QWidget *parent, bool withHeader) :
    QWidget(parent),
    ui(new Ui::MaterialPropertyTableWidget)
{
    ui->setupUi(this);

    if (withHeader) this->createHeader();
}

MaterialPropertyTableWidget::~MaterialPropertyTableWidget()
{
    delete ui;
}

void MaterialPropertyTableWidget::createHeader()
{
    QLabel* lblCol1 = new QLabel(tr("Name"));
    QLabel* lblCol2 = new QLabel(tr("Type"));
    QLabel* lblCol3 = new QLabel(tr("Definition"));
    QLabel* lblCol4 = new QLabel(tr("Unit"));
    QLabel* lblCol5 = new QLabel(tr("Symbol"));

    ui->grid->addWidget(lblCol1, 0, 0);
    ui->grid->addWidget(lblCol2, 0, 1);
    ui->grid->addWidget(lblCol3, 0, 2);
    ui->grid->addWidget(lblCol4, 0, 3);
    ui->grid->addWidget(lblCol5, 0, 4);
}

void MaterialPropertyTableWidget::addProperty(ECFParameter* property)
{
    this->addRow(
                QString::fromStdString(property->metadata.description.at(0)),
                property,
                QString::fromStdString(property->metadata.unit),
                QString::fromStdString(property->name)
                );
}

void MaterialPropertyTableWidget::addRow(const QString& name, ECFParameter* data,
                                         const QString& unit, const QString& symbol)
{
    int row = ui->grid->rowCount();

    QLabel* lblName = new QLabel(name, this);
    QLabel* lblUnit = new QLabel(unit, this);
    QLabel* lblSymbol = new QLabel(symbol, this);

    DataTypeEditWidget* dataWidget = new DataTypeEditWidget(data, this);
    QComboBox* cmbBox = dataWidget->createComboBox(this);
    //TODO: Variables

    ui->grid->addWidget(lblName, row, 0, Qt::AlignLeft| Qt::AlignTop);
    ui->grid->addWidget(cmbBox, row, 1, Qt::AlignTop);
    ui->grid->addWidget(dataWidget, row, 2, Qt::AlignTop);
    ui->grid->addWidget(lblUnit, row, 3, Qt::AlignLeft| Qt::AlignTop);
    ui->grid->addWidget(lblSymbol, row, 4, Qt::AlignLeft| Qt::AlignTop);

    this->m_properties.append(data);
    this->m_rowWidgets.append(dataWidget);
}

void MaterialPropertyTableWidget::save()
{
    foreach (DataTypeEditWidget* w, m_rowWidgets) {
        w->save();
    }
}

bool MaterialPropertyTableWidget::isValid()
{
    int index = 0;
    foreach (DataTypeEditWidget* w, m_rowWidgets) {
        if (!w->isValid())
        {
            this->m_invalidRow = index;
            return false;
        }
        index++;
    }

    return true;
}

QString MaterialPropertyTableWidget::errorMessage()
{
    if (m_rowWidgets.size() > 0)
    {
        return QString::fromStdString(m_properties.at(m_invalidRow)->metadata.description.at(0))
                + QLatin1String(": ")
                + m_rowWidgets.at(m_invalidRow)->errorMessage()
                + QLatin1String(".");
    }
    else
    {
        return QLatin1String("");
    }
}
