
#ifndef EXPRESSIONEDIT_H
#define EXPRESSIONEDIT_H

#include "ivalidatableobject.h"

#include <QLineEdit>
#include <QFocusEvent>

#include <vector>
#include <string>

namespace espreso
{

    class ExpressionEdit : public QLineEdit, public IValidatableObject
    {
        Q_OBJECT

    public:
        ExpressionEdit(const QString& contents,
                       const std::vector<std::string>& variables,
                       QWidget* parent = nullptr);

        bool isValid() override;
        QString errorMessage() override;

    private slots:
        void onTextChanged(const QString& text);

    private:
        bool m_isValid;
        std::vector<std::string> m_variables;
    };

}

#endif // EXPRESSIONEDIT_H
