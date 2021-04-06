#include "maptablewidget.h"

using namespace espreso;

MapTableWidget::MapTableWidget(const QStringList& headlines,
                               QWidget *parent,
                               std::unique_ptr<CleanRowFactory> rowFactory) :
    TableWidget(2, headlines, parent, std::move(rowFactory))
{

}

void MapTableWidget::addData(const QString& data)
{
    qWarning("MapTableWidget::addData: Empty method.");
}

QString MapTableWidget::data()
{
    return QLatin1String("");
}
