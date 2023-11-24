
#ifndef ECFOBJECTTREEWIDGET_H
#define ECFOBJECTTREEWIDGET_H

#include "config/configuration.h"

#include <QWidget>
#include <QStandardItem>
#include <QTreeView>
#include <QAction>

namespace espreso
{

struct ECFObjectTreeWidget : public QWidget
{
    Q_OBJECT

public:
    ECFObjectTreeWidget(const QString& label = "", QWidget* parent = 0);
    virtual ~ECFObjectTreeWidget() {}

    void add(ECFObject* obj);

protected:
    QTreeView* m_view;
    QStandardItemModel* m_model;
    QVector<QStandardItem*> m_groups;
    QVector<ECFObject*> m_objs;
    QStandardItem* m_root;

    virtual QDialog* createDialog(const QModelIndex& groupIndex, ECFParameter* param = nullptr) = 0;
    virtual QString dialogResult(QDialog* dialog) = 0;

	QModelIndex selectedGroupItem();
    ECFParameter* selectedParam(const QModelIndex& groupIndex);
    virtual std::string itemKeyInECFObject(QString nameInTree);

    virtual void newItemAccepted(int group, QString name) = 0;
    virtual void newItemRejected(int group) = 0;
	virtual void editItemAccepted(const QModelIndex& group, const QModelIndex& item, ECFParameter* param) = 0;
	virtual void editItemRejected(const QModelIndex& group, const QModelIndex& item, ECFParameter* param) = 0;
	virtual void deleteItemAccepted(const QModelIndex& group, int index, const QString& name) = 0;

private slots:
    void onActionNew();
    void onActionEdit();
    void onActionDelete();

    void onContextMenu(const QPoint &pos);

private:
    QAction* m_action_new;
    QAction* m_action_edit;
    QAction* m_action_delete;
};

}

#endif // ECFOBJECTTREEWIDGET_H
