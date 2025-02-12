
#ifndef FORMWIDGET_H
#define FORMWIDGET_H

#include "isavableobject.h"
#include "ivalidatableobject.h"
#include "fieldhandler.h"

#include "config/configuration.h"

#include <QWidget>
#include <QFormLayout>
#include <QLineEdit>
#include <QLabel>

namespace espreso
{

    class FormWidget : public QWidget, public ISavableObject,
            public IValidatableObject
    {
        Q_OBJECT

    public:
        explicit FormWidget(QWidget* = nullptr);

        void appendString(ECFParameter*);
        void appendNonnegativeInteger(ECFParameter*);
        void appendPositiveInteger(ECFParameter*);
        void appendFloat(ECFParameter*);

        void appendRow(const QString&, QWidget*);

        void save() override;
        void saveState() override;
        void restoreState() override;

        bool isValid() override;
        QString errorMessage() override;

    private:
        QFormLayout* m_layout;

        QVector<QPair<ECFParameter*, FieldHandler*> > m_fields;
        int m_invalidIndex = 0;

        QVector<std::string> m_state_fields;
        bool m_stateStored = false;

        QLabel* createLabel(const QString& text);
        void appendRow(QLabel*, QWidget*);
    };

}

#endif // FORMWIDGET_H
