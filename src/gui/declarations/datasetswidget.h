
#ifndef DATASETSWIDGET_H
#define DATASETSWIDGET_H

#include "material/materialdialog.h"
#include "gui/elements/ecfobjecttreewidget.h"

namespace espreso
{

class DataSetsWidget : public ECFObjectTreeWidget
{
    Q_OBJECT
public:
    DataSetsWidget(ECFObject* materials, const QString& label = "", QWidget* parent = 0);

    void setMaterials(ECFObject* materials);

protected:
    virtual QDialog* createDialog(const QModelIndex& groupIndex, ECFParameter* param = nullptr) override;
    virtual QString dialogResult(QDialog* dialog) override;
	virtual void newItemAccepted(int, QString) override {}
    virtual void newItemRejected(int group) override;
	virtual void editItemAccepted(const QModelIndex&, const QModelIndex&, ECFParameter*) override {}
    virtual void editItemRejected(const QModelIndex &, const QModelIndex &, ECFParameter *) override {}
	virtual void deleteItemAccepted(const QModelIndex& group, int index, const QString& name) override;
    virtual std::string itemKeyInECFObject(QString nameInTree) override;

private:
    ECFObject* m_materials;
    int m_materials_id = 1;
    QVector<std::string> m_materials_names;
    QVector<std::string> m_materials_ids;
    MaterialConfiguration* m_last_modified;

    void initMaterials();
    MaterialConfiguration* newMaterial();
};

}

#endif // DATASETSWIDGET_H
