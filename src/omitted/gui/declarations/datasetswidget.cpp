#include "datasetswidget.h"

#include <QDebug>

using namespace espreso;

DataSetsWidget::DataSetsWidget(ECFObject* materials,
                               const QString& label,
                               QWidget* parent) :
    ECFObjectTreeWidget(label, parent)
{
    this->m_materials = materials;
    this->initMaterials();
}

void DataSetsWidget::initMaterials()
{
    QStandardItem* material_group = new QStandardItem(QString::fromStdString(this->m_materials->name));
    this->m_root->appendRow(material_group);
    this->m_groups.append(material_group);
    this->m_objs.append(this->m_materials);
    this->m_view->setExpanded(material_group->index(), true);

    for (auto m = this->m_materials->parameters.begin();
         this->m_materials->parameters.end() != m;
         ++m)
    {
        bool ok;
        int mid = QString::fromStdString((*m)->name).toInt(&ok);
        if (ok && mid > this->m_materials_id) this->m_materials_id = mid;

        ECFObject* material = static_cast<ECFObject*>(*m);

        this->m_materials_names.append(material->getParameter("name")->getValue());
        this->m_materials_ids.append((*m)->name);

        QStandardItem* item = new QStandardItem(
                    QString::fromStdString(
                        this->m_materials_names.last()
                        )
                    );
        material_group->appendRow(item);
    }

    this->m_materials_id++;
}

void DataSetsWidget::setMaterials(ECFObject* materials)
{
    if (this->m_materials == materials) return;

    materials->dropAllParameters();

    for (auto m = this->m_materials->parameters.begin();
         m != this->m_materials->parameters.end();
         ++m)
    {
        materials->parameters.push_back(*m);
    }

    this->m_materials = materials;
}

QDialog* DataSetsWidget::createDialog(const QModelIndex& groupIndex, ECFParameter* param)
{
    //ECFObject* obj = this->m_objs[groupIndex.row()];

    if (groupIndex.row() == 0)
    {
        MaterialConfiguration* mc;

        if (param == nullptr) mc = this->newMaterial();
        else mc = static_cast<MaterialConfiguration*>(param);

        int index = this->m_materials_names.indexOf(mc->getParameter("name")->getValue());
        if (index >= 0)
        {
            this->m_materials_names.remove(index);
            this->m_materials_ids.remove(index);
        }

        MaterialDialog* md = new MaterialDialog(mc, this->m_materials_names, this);
        this->m_last_modified = mc;

        return md;
    }

    qFatal("DataSetsWidget: Failed to create dialog. Unknown object group: %d", groupIndex.row());
    return nullptr;
}

QString DataSetsWidget::dialogResult(QDialog*)
{
    std::string name = this->m_last_modified->name;
	this->m_materials_names.append(name);
	this->m_materials_ids.append(std::to_string(this->m_materials_id - 1));
    return QString::fromStdString(name);
}

void DataSetsWidget::deleteItemAccepted(const QModelIndex& group, int, const QString& name)
{
	if (group.row() == 0)
	{
		int matIndex = this->m_materials_names.indexOf(name.toStdString());
		this->m_materials_names.remove(matIndex);
		this->m_materials_ids.remove(matIndex);
	}
}

MaterialConfiguration* DataSetsWidget::newMaterial()
{
    return static_cast<MaterialConfiguration*>(
                this->m_materials->getParameter(std::to_string(this->m_materials_id++))
                );
}

void DataSetsWidget::newItemRejected(int group)
{
    if (group == 0)
    {
        MaterialConfiguration* tmp = this->m_last_modified;
        this->m_last_modified = nullptr;
        this->m_materials_id--;

        this->m_materials->dropParameter(tmp);
    }
}

std::string DataSetsWidget::itemKeyInECFObject(QString nameInTree)
{
    int index = this->m_materials_names.indexOf(nameInTree.toStdString());

    return this->m_materials_ids[index];
}
