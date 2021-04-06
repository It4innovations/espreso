#include "expressioneditdelegate.h"

#include "basis/expression/expression.h"

using namespace espreso;

ExpressionEditDelegate::ExpressionEditDelegate(const std::vector<std::string>& variables,
                                               QObject *parent)
    : QItemDelegate(parent)
{
    this->m_variables = variables;
}

QWidget* ExpressionEditDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QLineEdit *editor = new QLineEdit(parent);

    return editor;
}


void ExpressionEditDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    line->setText(value);
}


void ExpressionEditDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index, value);
}


void ExpressionEditDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

void ExpressionEditDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    QVariant data = index.data();
    QString text = data.toString();
    QStyleOptionViewItem newOption(option);
    if (!Expression::isValid(text.toStdString(), m_variables))
    {
        newOption.font.setBold(true);
        newOption.palette.setColor(QPalette::Text, Qt::red);
    }

    QItemDelegate::paint(painter, newOption, index);
}
