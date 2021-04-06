#ifndef ECFPARAMETERTREEWIDGET_H
#define ECFPARAMETERTREEWIDGET_H

#include <QWidget>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QItemDelegate>
#include <functional>

#include "config/configuration.h"
#include "isavableobject.h"
#include "ivalidatableobject.h"
#include "ecfparametertreedelegate.h"

namespace espreso {

namespace Ui {
class ECFParameterTreeWidget;
}

class ECFParameterTreeWidget : public QWidget,
        public IValidatableObject,
        public ISavableObject
{
    Q_OBJECT

public:
    explicit ECFParameterTreeWidget(QWidget *parent = 0);
    ~ECFParameterTreeWidget();

    int createGroup(const QString& name, int parent_id = 0);

    void add(ECFValue *value, int parent_id = 0);
    void addString(ECFValue *string, int parent_id = 0);
    void addNonnegativeInteger(ECFValue* nnint, int parent_id = 0);
    void addPositiveInteger(ECFValue* pint, int parent_id = 0);
    void addFloat(ECFValue* p_float, int parent_id = 0);
    void addOption(ECFValue* option, int parent_id = 0);
    void addBool(ECFValue* p_bool, int parent_id = 0);
    void addWithEditorFactory(ECFValue* value, TextItemWidgetFactory* editorFactory, int parent_id = 0);
    void addWithWidget(ECFValue* value, QWidget *widget, int parent_id = 0);
    void resizeCellsToContent();

    virtual bool isValid() override;
    virtual QString errorMessage() override;
    virtual void save() override;
    virtual void saveState() override;
    virtual void restoreState() override;

signals:
    void itemChanged();

private:
    Ui::ECFParameterTreeWidget *ui;

    QStandardItemModel* m_model;
    QVector<ECFValue*> m_values;
    QStandardItem* m_root;
    QVector<QStandardItem*> m_groups;

    ECFParameterTreeDelegate* m_delegate;

    QVector<std::string> m_stored_values;

    QString m_err;

    void addRow(ECFValue* value, TextItemWidgetFactory* editorFactory, int parent_id, bool valueCheckable = false);
    void applyForAllItems(QStandardItem* root, std::function<bool(QStandardItem*)> apply);
    bool findIndex(QStandardItem* item, QStandardItem* root, int *index);

private slots:
    void onItemChanged(QStandardItem* item);
};

}

#endif // ECFPARAMETERTREEWIDGET_H
