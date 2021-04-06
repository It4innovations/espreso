#ifndef REGIONPROPERTYWIDGET_H
#define REGIONPROPERTYWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QTreeView>

#include "regionpairdialog.h"
#include "ecfobjecttreewidget.h"

namespace espreso
{

struct ECFObject;
struct ECFParameter;
class Mesh;

class RegionPropertyWidget : public ECFObjectTreeWidget
{
    Q_OBJECT

public:
    explicit RegionPropertyWidget(Mesh* mesh, PhysicsConfiguration* physics, const QString& label = "", QWidget *parent = 0);

    void addProperty(ECFObject* obj);

protected:
    virtual QDialog* createDialog(const QModelIndex& groupIndex, ECFParameter* param = nullptr) override;
    virtual QString dialogResult(QDialog* dialog) override;

    virtual void newItemAccepted(int, QString) override {}
    virtual void newItemRejected(int) override {}
	virtual void editItemAccepted(const QModelIndex&, const QModelIndex&, ECFParameter*) override {}
	virtual void editItemRejected(const QModelIndex &group, const QModelIndex &item, ECFParameter *param) override {}
	virtual void deleteItemAccepted(const QModelIndex& group, int index, const QString& name) override {}

private:
    Mesh* m_mesh;
    PhysicsConfiguration* m_physics;
};

}

#endif // REGIONPROPERTYWIDGET_H
