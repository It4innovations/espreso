
#include "regionpropertywidget.h"

#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"

#include <QMenu>
#include <QTreeView>
#include <QDebug>

using namespace espreso;

RegionPropertyWidget::RegionPropertyWidget(const QString& label, QWidget *parent) :
    ECFObjectTreeWidget(label, parent)
{

}

void RegionPropertyWidget::addProperty(ECFObject *obj)
{
    this->add(obj);
}

QDialog* RegionPropertyWidget::createDialog(const QModelIndex& groupIndex, ECFParameter *param)
{
    if (m_objs[groupIndex.row()]->metadata.datatype.size() == 2)
    {
        if (m_objs[groupIndex.row()]->metadata.datatype[1] == ECFDataType::EXPRESSION)
            return RegionPairDialog::createRegionExpression(m_objs[groupIndex.row()], info::mesh, param);

        if (m_objs[groupIndex.row()]->metadata.datatype[1] == ECFDataType::MATERIAL)
            return RegionPairDialog::createRegionMaterial(m_objs[groupIndex.row()], info::mesh,
                    static_cast<ECFObject*>(info::ecf->getPhysics()->ecfdescription->getParameter(&info::ecf->getPhysics()->materials)), param);

        if (m_objs[groupIndex.row()]->metadata.datatype[1] == ECFDataType::OPTION)
            return RegionPairDialog::createRegionOption(m_objs[groupIndex.row()], info::mesh, param);
    }
    if (m_objs[groupIndex.row()]->metadata.datatype.size() == 1
            && m_objs[groupIndex.row()]->metadata.description.size() == 2)
    {
        ECFObject* pair = nullptr;
        if (param != nullptr)
        {
            if (!param->isObject())
            {
                qCritical("RegionPairDialog: Failed to create a dialog. ECFObject* expected. ECFValue* given.");
                return nullptr;
            }
            else {
                pair = static_cast<ECFObject*>(param);
            }
        }
        return RegionPairDialog::createRegionObject(m_objs[groupIndex.row()], info::mesh, pair);
    }

    qCritical("RegionPairDialog: Failed to create a dialog. Unknown type of selected object.");
    return nullptr;
}

QString RegionPropertyWidget::dialogResult(QDialog* dialog)
{
    return static_cast<RegionPairDialog*>(dialog)->region();
}
