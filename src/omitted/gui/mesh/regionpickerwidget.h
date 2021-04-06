#ifndef REGIONPICKERWIDGET_H
#define REGIONPICKERWIDGET_H

#include <QWidget>
#include "meshwidget.h"

namespace espreso
{

    namespace Ui {
    class RegionPickerWidget;
    }

    class RegionPickerWidget : public QWidget
    {
        Q_OBJECT

    public:
        explicit RegionPickerWidget(MeshWidget* mesh, QWidget *parent = 0);
        ~RegionPickerWidget();

    private:
        Ui::RegionPickerWidget *ui;

        QStandardItemModel* m_model;
        void setupModel();

        MeshWidget* m_mesh;

        QVector<QString> m_regions;

    private slots:
        void regionStateChanged(QStandardItem* item);
        void onRegionClicked(const QString& region);
    };

}

#endif // REGIONPICKERWIDGET_H
