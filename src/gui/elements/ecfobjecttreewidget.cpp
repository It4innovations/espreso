
#include "ecfobjecttreewidget.h"

#include <QVBoxLayout>
#include <QLabel>
#include <QDialog>
#include <QMenu>
#include <QDebug>

using namespace espreso;

ECFObjectTreeWidget::ECFObjectTreeWidget(const QString& label, QWidget* parent) :
    QWidget(parent)
{
    QVBoxLayout* layout = new QVBoxLayout;
    layout->setMargin(0);
    layout->setSpacing(0);
    this->setLayout(layout);

    if (!label.isEmpty())
    {
        QLabel* headline = new QLabel(label, this);
        layout->addWidget(headline);
    }

    this->m_model = new QStandardItemModel();
    this->m_root = this->m_model->invisibleRootItem();

    QTreeView* view = new QTreeView(this);
    view->setModel(m_model);
    view->setEditTriggers(QAbstractItemView::NoEditTriggers);
    view->setContextMenuPolicy(Qt::CustomContextMenu);
    view->setHeaderHidden(true);
    layout->addWidget(view);
    connect(view, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(onContextMenu(QPoint)));
    this->m_view = view;

    this->m_action_new = new QAction(tr("&New"), this);
    connect(this->m_action_new, SIGNAL(triggered()), this, SLOT(onActionNew()));

    this->m_action_edit = new QAction(tr("&Edit"), this);
    connect(this->m_action_edit, SIGNAL(triggered()), this, SLOT(onActionEdit()));

    this->m_action_delete = new QAction(tr("&Delete"), this);
    connect(this->m_action_delete, SIGNAL(triggered()), this, SLOT(onActionDelete()));
}

void ECFObjectTreeWidget::add(ECFObject *obj)
{
    QStandardItem* group = new QStandardItem(QString::fromStdString(obj->name));
    this->m_root->appendRow(group);
    this->m_groups.append(group);
    this->m_objs.append(obj);

    this->m_view->setExpanded(group->index(), true);

    for (auto it = obj->parameters.begin(); it != obj->parameters.end(); ++it)
    {
        if ((*it)->name.compare("INTERSECTION") == 0) continue;
        group->appendRow(new QStandardItem(QString::fromStdString((*it)->name)));
    }
}

void ECFObjectTreeWidget::onActionNew()
{
    QModelIndex groupIndex = this->selectedGroupItem();
    if ( !(groupIndex.isValid()) ) return;

    QDialog* dialog = this->createDialog(groupIndex);

    if (dialog->exec() == QDialog::Accepted)
    {
        QString item_name = this->dialogResult(dialog);
        QStandardItem* item = new QStandardItem(item_name);
        this->m_groups[groupIndex.row()]->appendRow(item);
        this->newItemAccepted(groupIndex.row(), item_name);
    }
    else
    {
        this->newItemRejected(groupIndex.row());
    }
}

void ECFObjectTreeWidget::onActionEdit()
{
    QModelIndex groupIndex = this->selectedGroupItem();
    if ( !(groupIndex.isValid()) ) return;

    ECFParameter* param = this->selectedParam(groupIndex);
    if (param == nullptr) return;

    QModelIndexList indexList = m_view->selectionModel()->selectedIndexes();
    QModelIndex clicked = indexList.at(0);

    QDialog* dialog = this->createDialog(groupIndex, param);

    int res = dialog->exec();
    QString item_name = this->dialogResult(dialog);
    this->m_model->setData(clicked, item_name);

    if (res == QDialog::Accepted) this->editItemAccepted(groupIndex, clicked, param);
    else this->editItemRejected(groupIndex, clicked, param);
}

void ECFObjectTreeWidget::onActionDelete()
{
    QModelIndex groupIndex = this->selectedGroupItem();
    if ( !(groupIndex.isValid()) ) return;

    ECFParameter* param = this->selectedParam(groupIndex);
    if (param == nullptr) return;

    this->m_objs[groupIndex.row()]->dropParameter(param);

    QModelIndexList indexList = m_view->selectionModel()->selectedIndexes();
    QModelIndex clicked = indexList.at(0);

    int itemIndex = clicked.row();
    QString itemName = clicked.data().toString();

    this->m_groups[groupIndex.row()]->removeRow(clicked.row());

    this->deleteItemAccepted(groupIndex, itemIndex, itemName);
}

ECFParameter* ECFObjectTreeWidget::selectedParam(const QModelIndex &groupIndex)
{
    QModelIndexList indexList = m_view->selectionModel()->selectedIndexes();
    QModelIndex clicked = indexList.at(0);
    if (!clicked.parent().isValid()) return nullptr;

    ECFObject* obj = this->m_objs[groupIndex.row()];
    QString item_name =  this->m_model->data(clicked).toString();
    std::string key = this->itemKeyInECFObject(item_name);

    return obj->getParameter(key);
}

std::string ECFObjectTreeWidget::itemKeyInECFObject(QString nameInTree)
{
    return nameInTree.toStdString();
}

QModelIndex ECFObjectTreeWidget::selectedGroupItem()
{
    QModelIndexList indexList = m_view->selectionModel()->selectedIndexes();
    if (!indexList.count())
        return QModelIndex();

    QModelIndex clicked = indexList.at(0);
    QModelIndex parent = clicked.parent();

    if (!parent.isValid())
    {
        parent = clicked;
    }

    return parent;
}

void ECFObjectTreeWidget::onContextMenu(const QPoint &pos)
{
    QMenu treeMenu(this);
    treeMenu.addAction(this->m_action_new);
    treeMenu.addAction(this->m_action_edit);
    treeMenu.addAction(this->m_action_delete);
    treeMenu.exec(m_view->mapToGlobal(pos));
}
