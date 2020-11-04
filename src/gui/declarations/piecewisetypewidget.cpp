
#include "piecewisetypewidget.h"

#include "gui/elements/expressioneditdelegate.h"
#include "gui/elements/expressionedit.h"
#include "gui/elements/comboboxdelegate.h"
#include "gui/validators/validatordelegate.h"

#include "basis/expression/expression.h"

#include <QDebug>

using namespace espreso;

QStringList PiecewiseTypeWidget::headlines()
{
    QStringList result;
    result << QString("O/C")
           << QObject::tr("Lower bound")
           << QObject::tr("Upper bound")
           << QString("O/C")
           << QObject::tr("f(x)");

    return result;
}

PiecewiseTypeWidget::PiecewiseTypeWidget(const std::vector<std::string>& variables,
                                         QWidget* parent) :
    TableWidget(5, PiecewiseTypeWidget::headlines(), parent)
{
    this->m_variables = variables;
    QStringList leftCmbOptions;
    leftCmbOptions << "<" << "(";
    QStringList rightCmbOptions;
    rightCmbOptions << ">" << ")";
    this->mTable->setItemDelegateForColumn(0, new ComboBoxDelegate(leftCmbOptions, this));
    this->mTable->setItemDelegateForColumn(1, new ValidatorDelegate(new DoubleValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(2, new ValidatorDelegate(new DoubleValidatorFactory(), this));
    this->mTable->setItemDelegateForColumn(3, new ComboBoxDelegate(rightCmbOptions, this));
    this->mTable->setItemDelegateForColumn(4, new ExpressionEditDelegate(variables, this));
}

bool PiecewiseTypeWidget::isValid()
{
    for (int row = 0; row < mModel->rowCount() - 1; ++row)
    {
        QModelIndex index = this->mModel->index(row, 4);
        QString expression = index.data().toString();
        if (!Expression::isValid(expression.toStdString(), m_variables))
        {
            this->m_invalidRow = row + 1;
            this->m_isValid = false;
            return false;
        }
    }

    this->m_isValid = true;

    return TableWidget::isValid();
}

QString PiecewiseTypeWidget::errorMessage()
{
    if (m_isValid)
        return TableWidget::errorMessage();
    else
        return tr("Invalid f(x) formula at row %1")
                .arg(m_invalidRow);
}


void PiecewiseTypeWidget::addData(const QString& data)
{
    QRegExp rxCase("IF");
    QStringList ifs = data.split(rxCase, QString::SkipEmptyParts);

    foreach (QString v_if, ifs) {

        QVector<QString> row;

        QStringList middle = v_if.split(')', QString::SkipEmptyParts);

        QString condition = middle.at(0);
        condition.remove(0, 1);
        QStringList bounds = condition.split("AND", QString::SkipEmptyParts);

        QString leftBound = bounds.at(0);
        QString rightBound = bounds.at(1);

        QString left;
        if (leftBound.contains(">="))
        {
            left = leftBound.split(">=", QString::SkipEmptyParts).at(1);
            row.append("<");
        }
        else if (leftBound.contains('>'))
        {
            left = leftBound.split(">", QString::SkipEmptyParts).at(1);
            row.append("(");
        }
        row.append(left.trimmed());

        QString right;
        if (rightBound.contains("<="))
        {
            right = rightBound.split("<=", QString::SkipEmptyParts).at(1);
            row.append(right.trimmed());
            row.append(">");
        }
        else if (rightBound.contains('<'))
        {
            right = rightBound.split('<', QString::SkipEmptyParts).at(1);
            row.append(right.trimmed());
            row.append(")");
        }

        middle.removeFirst();
        QString expr = middle.join(")");
        expr = expr.split(';', QString::SkipEmptyParts).at(0);
        expr = expr.trimmed();
        row.append(expr);

        this->addRow(row);
    }
}

QString PiecewiseTypeWidget::data()
{
    QStringList ifs;
    for (int i = 0; i < mModel->rowCount() - 1; ++i)
    {
        QLatin1String start("IF (X");
        QModelIndex left = mModel->index(i, 0);
        QModelIndex lBound = mModel->index(i, 1);
        QString leftSign;
        if (mModel->data(left).toString().compare(QLatin1String("(")) == 0)
            leftSign = QLatin1String(" > ");
        else
            leftSign = QLatin1String(" >= ");
        QModelIndex rBound = mModel->index(i, 2);
        QModelIndex right = mModel->index(i, 3);
        QString rightSign;
        if (mModel->data(right).toString().compare(QLatin1String(")")) == 0)
            rightSign = QLatin1String(" < ");
        else
            rightSign = QLatin1String(" <= ");
        QModelIndex fnIndex = mModel->index(i, 4);
        QString line = start
                + leftSign + mModel->data(lBound).toString()
                + QLatin1String(" AND X")
                + rightSign + mModel->data(rBound).toString()
                + QLatin1String(") ") + mModel->data(fnIndex).toString()
                + QLatin1String("; ");
        ifs << line;
    }

    return ifs.join(QLatin1String(" "));
}
