#ifndef TABLEWIDGET_H
#define TABLEWIDGET_H

#include "gui/elements/ivalidatableobject.h"
#include "gui/validators/validatorfactory.h"

#include <QWidget>
#include <QStandardItemModel>
#include <QStringList>
#include <QTableView>
#include <QAction>
#include <QMenu>

#include <memory>

namespace espreso
{

    namespace Ui {
    class TableWidget;
    }

    class CleanRowFactory
    {
    public:
        virtual QList<QStandardItem*> create(int columns) = 0;
    };

    class DefaultCleanRowFactory : public CleanRowFactory
    {
    public:
        virtual QList<QStandardItem*> create(int columns) override;
    };

    class TableWidget : public QWidget, public IValidatableObject
    {
        Q_OBJECT

    public:
        virtual ~TableWidget();

        virtual void addRow(const QVector<QString>& rowData);
        virtual void addData(const QVector<QVector<QString> >& data);
        virtual void addData(const QString& data) = 0;
        virtual bool isValid() override;
        virtual QString errorMessage() override;
        virtual QString data() = 0;

    signals:
        void cellChanged();

    protected:
        QTableView* mTable;
        QStandardItemModel* mModel;

        explicit TableWidget(int columns,
                             const QStringList& headlines,
                             QWidget *parent = 0,
                             std::unique_ptr<CleanRowFactory> rowFactory = std::unique_ptr<CleanRowFactory>(new DefaultCleanRowFactory));

    private slots:
        void onItemDoubleClick(const QModelIndex& index);
        void onContextMenu(const QPoint& point);
        void deleteItem();
        void onDataChanged(QModelIndex, QModelIndex, QVector<int>);

    private:
        Ui::TableWidget *ui;

        std::unique_ptr<CleanRowFactory> m_cleanRowFactory;

        QModelIndex toDelete;
        QAction* mActionDelete;

        int m_invalidRow = 0;
        int m_invalidCol = 0;
        bool m_empty = true;
    };
}
#endif // TABLEWIDGET_H
