#include "materialwidget.h"

#include <QVBoxLayout>
#include <QScrollArea>
#include <QComboBox>
#include <QMessageBox>
#include <QLabel>
#include <QScrollArea>
#include <QLineEdit>
#include <QDebug>
#include "elements/optionhandler.h"
#include "elements/formwidget.h"
#include "elements/boolhandler.h"
#include "materialpropertytablewidget.h"
#include "utils/utils.h"

using namespace espreso;

MaterialWidget::MaterialWidget(MaterialConfiguration* material,
                               const QVector<std::string>& materialNames,
                               QWidget *parent) :
    ScrollECFObjectWidget(material, parent)
{
    this->m_names = materialNames;
}

bool MaterialWidget::isValid()
{
    if (!ECFObjectWidget::isValid())
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(this->errorMessage());
        msg.exec();

        return false;
    }

    ISavableObject* info = m_savables.first();
    info->saveState();
    info->save();

    std::string newName =
            QString::fromStdString(
                m_obj->getParameter(std::string("name"))->getValue()
                )
                .toUpper()
                .toStdString();

    if (Utils::caseInsensitiveIndexOf(m_names, newName) != -1)
    {
        QMessageBox msg;
        msg.setWindowTitle(tr("Error"));
        msg.setText(tr("Material with same name already exists!"));
        msg.exec();

        info->restoreState();

        return false;
    }

    foreach (ISavableObject* obj, m_savables) {
        obj->save();
    }

    return true;
}
