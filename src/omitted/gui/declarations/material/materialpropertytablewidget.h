#ifndef MATERIALPROPERTYTABLEWIDGET_H
#define MATERIALPROPERTYTABLEWIDGET_H

#include <QWidget>
#include <QVector>
#include "config/configuration.h"
#include "elements/isavableobject.h"
#include "elements/ivalidatableobject.h"
#include "datatypeeditwidget.h"

namespace espreso
{

    namespace Ui {
    class MaterialPropertyTableWidget;
    }

    class MaterialPropertyTableWidget : public QWidget,
            public ISavableObject, public IValidatableObject
    {
        Q_OBJECT

    public:
        explicit MaterialPropertyTableWidget(QWidget *parent = 0,
                                             bool withHeader = true);
        ~MaterialPropertyTableWidget();

        void addProperty(ECFParameter* property);
        void addRow(const QString& name, ECFParameter* data,
                    const QString& unit, const QString& symbol);

        void save() override;

        bool isValid() override;
        QString errorMessage() override;

    private:
        Ui::MaterialPropertyTableWidget *ui;

        QVector<ECFParameter*> m_properties;
        QVector<DataTypeEditWidget*> m_rowWidgets;

        int m_invalidRow = 0;

        void createHeader();
    };

}
#endif // MATERIALPROPERTYTABLEWIDGET_H
