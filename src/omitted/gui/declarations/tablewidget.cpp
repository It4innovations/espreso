#include "tablewidget.h"
#include "ui_tablewidget.h"

#include <QtDebug>

using namespace espreso;

QList<QStandardItem*> DefaultCleanRowFactory::create(int columns)
{
    QList<QStandardItem*> row;
    for (int i = 0; i < columns; ++i)
    {
        row << new QStandardItem();
    }

    return row;
}

TableWidget::TableWidget(int columns, const QStringList& headlines,
                         QWidget *parent, std::unique_ptr<CleanRowFactory> rowFactory) :
    QWidget(parent),
    ui(new Ui::TableWidget)
{
    ui->setupUi(this);

    this->mTable = this->ui->table;

    this->mModel = new QStandardItemModel(this);
    ui->table->setModel(mModel);

    this->mModel->setColumnCount(columns);

    this->mModel->setHorizontalHeaderLabels(headlines);

    this->m_cleanRowFactory = std::move(rowFactory);
    this->mModel->appendRow(this->m_cleanRowFactory->create(this->mModel->columnCount()));

    connect(mTable, SIGNAL(doubleClicked(QModelIndex)),
            this, SLOT(onItemDoubleClick(QModelIndex)));

    this->mTable->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(mTable, SIGNAL(customContextMenuRequested(QPoint)),
            this, SLOT(onContextMenu(QPoint)));

    this->mActionDelete = new QAction(tr("&Delete"), this);
    connect(mActionDelete, SIGNAL(triggered(bool)), this, SLOT(deleteItem()));

	connect(this->mModel, SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),
								 this, SLOT(onDataChanged(QModelIndex,QModelIndex,QVector<int>)));
}

TableWidget::~TableWidget()
{
    delete ui;
}

void TableWidget::addRow(const QVector<QString>& rowData)
{
    QList<QStandardItem*> row;
    foreach (QString data, rowData) {
        row.append(new QStandardItem(data));
    }
    this->mModel->insertRow(mModel->rowCount() - 1, row);
}

void TableWidget::addData(const QVector<QVector<QString> >& data)
{
    foreach (QVector<QString> row, data) {
        this->addRow(row);
    }
}

void TableWidget::onItemDoubleClick(const QModelIndex &index)
{
    if (index.row() != this->mModel->rowCount() - 1)
        return;

    this->mModel->appendRow(this->m_cleanRowFactory->create(this->mModel->columnCount()));
}

void TableWidget::onContextMenu(const QPoint& point)
{
    this->toDelete = this->mTable->indexAt(point);
    QMenu menu(this);
    menu.addAction(this->mActionDelete);
    menu.exec(ui->table->mapToGlobal(point));
}

void TableWidget::deleteItem()
{
    if (toDelete.row() == -1 || toDelete.column() == -1)
        return;
    if (toDelete.row() == mModel->rowCount() - 1)
        return;

    this->mModel->removeRow(toDelete.row());
}

bool TableWidget::isValid()
{
    if (mModel->rowCount() <= 1)
    {
        this->m_empty = true;
        return false;
    }

    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        for (int col = 0; col < mModel->columnCount(); ++col)
        {
            QModelIndex index = this->mModel->index(row, col);
            QString val = index.data().toString();
            if (val.isEmpty())
            {
                this->m_invalidRow = row + 1;
                this->m_invalidCol = col + 1;
                this->m_empty = false;
                return false;
            }
        }
    }

    return true;
}

QString TableWidget::errorMessage()
{
    if (m_empty)
        return tr("Empty table");
    else
        return tr("Empty cell at row %1 and column %2")
                .arg(m_invalidRow)
                .arg(m_invalidCol);
}

void TableWidget::onDataChanged(QModelIndex, QModelIndex, QVector<int>)
{
	emit cellChanged();
}
