
#include "regionpickerwidget.h"
#include "ui_regionpickerwidget.h"

using namespace espreso;

RegionPickerWidget::RegionPickerWidget(MeshWidget* mesh, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::RegionPickerWidget)
{
    ui->setupUi(this);

    this->m_mesh = mesh;
    connect(m_mesh, SIGNAL(regionClicked(QString)), this, SLOT(onRegionClicked(QString)));

    this->setupModel();
}

RegionPickerWidget::~RegionPickerWidget()
{
    delete ui;
}

void RegionPickerWidget::setupModel()
{
    this->m_model = new QStandardItemModel(ui->tree);

    QList<QString> regions = m_mesh->regions();
    m_regions = QVector<QString>::fromList(regions);

    foreach (QString r, regions) {
        QStandardItem* item = new QStandardItem(r);
        item->setCheckable(true);
        item->setCheckState(Qt::Checked);
        m_model->appendRow(item);
    }

    ui->tree->setModel(m_model);

    connect(m_model, SIGNAL(itemChanged(QStandardItem*)),
            this, SLOT(regionStateChanged(QStandardItem*)));
}

void RegionPickerWidget::regionStateChanged(QStandardItem* item)
{
    this->m_mesh->changeRegionState(item->text());
}

void RegionPickerWidget::onRegionClicked(const QString& region)
{
    int index = m_regions.indexOf(region);
    QModelIndex idx = m_model->index(index, 0);
    ui->tree->setCurrentIndex(idx);
}
