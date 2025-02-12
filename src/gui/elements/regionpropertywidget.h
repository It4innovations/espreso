#ifndef REGIONPROPERTYWIDGET_H
#define REGIONPROPERTYWIDGET_H

#include "regionpairdialog.h"
#include "ecfobjecttreewidget.h"

#include <QWidget>
#include <QStandardItemModel>
#include <QTreeView>

namespace espreso
{

struct ECFObject;
struct ECFParameter;

class RegionPropertyWidget : public ECFObjectTreeWidget
{
    Q_OBJECT

public:
    explicit RegionPropertyWidget(const QString& label = "", QWidget *parent = 0);

    void addProperty(ECFObject* obj);

protected:
    virtual QDialog* createDialog(const QModelIndex& groupIndex, ECFParameter* param = nullptr) override;
    virtual QString dialogResult(QDialog* dialog) override;

    virtual void newItemAccepted(int, QString) override {}
    virtual void newItemRejected(int) override {}
    virtual void editItemAccepted(const QModelIndex&, const QModelIndex&, ECFParameter*) override {}
    virtual void editItemRejected(const QModelIndex &group, const QModelIndex &item, ECFParameter *param) override {}
    virtual void deleteItemAccepted(const QModelIndex& group, int index, const QString& name) override {}
};

}

#endif // REGIONPROPERTYWIDGET_H
