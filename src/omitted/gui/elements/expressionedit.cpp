#include "expressionedit.h"

#include "basis/expression/expression.h"

using namespace espreso;


ExpressionEdit::ExpressionEdit(const QString& contents,
                               const std::vector<std::string>& variables,
                               QWidget* parent) :
    QLineEdit(contents, parent)
{
    this->m_variables = variables;

    connect(this, SIGNAL(textChanged(QString)),
            this, SLOT(onTextChanged(QString)));

    this->m_isValid = Expression::isValid(contents.toStdString(), m_variables);
}

bool ExpressionEdit::isValid()
{
    return this->m_isValid;
}

QString ExpressionEdit::errorMessage()
{
    if (this->text().isEmpty())
        return tr("Empty expression");
    return tr("Invalid expression");
}

void ExpressionEdit::onTextChanged(const QString& text)
{
    this->m_isValid = Expression::isValid(text.toStdString(), m_variables);
}
