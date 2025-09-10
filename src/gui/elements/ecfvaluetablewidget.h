
#ifndef ECFVALUETABLEWIDGET_H
#define ECFVALUETABLEWIDGET_H

#include "isavableobject.h"
#include "ivalidatableobject.h"

#include "config/configuration.h"

#include <QWidget>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QItemDelegate>

namespace espreso
{

namespace Ui {
struct ECFValueTableWidget;
}

struct ECFValueTableWidget : public QWidget,
        public IValidatableObject,
        public ISavableObject
{
    Q_OBJECT

public:
    explicit ECFValueTableWidget(QWidget *parent = 0);
    ~ECFValueTableWidget();

    void add(ECFValue *value);
    void addString(ECFValue *string);
    void addNonnegativeInteger(ECFValue* nnint);
    void addPositiveInteger(ECFValue* pint);
    void addFloat(ECFValue* p_float);
    void addOption(ECFValue* option);
    void addBool(ECFValue* p_bool);
    void addWithDelegate(ECFValue* value, QItemDelegate* delegate);
    void addWithWidget(ECFValue* value, QWidget *widget);
    void resizeCellsToContent();

    bool isValid() override;
    QString errorMessage() override;
    void save() override;
    void saveState() override;
    void restoreState() override;

signals:
    void itemChanged();

private:
    Ui::ECFValueTableWidget *ui;

    QStandardItemModel* m_model;
    QVector<ECFValue*> m_values;

    QVector<std::string> m_stored_values;

    int m_invalid = -1;

private slots:
    void onItemChanged(QStandardItem* item);
};

}
#endif // ECFVALUETABLEWIDGET_H
