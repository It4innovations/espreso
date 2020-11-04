
#include "validatordelegate.h"

#include <QLineEdit>

ValidatorDelegate::ValidatorDelegate(ValidatorFactory* factory, QObject *parent) : QItemDelegate(parent)
{
    this->mFactory = factory;
}
ValidatorDelegate::~ValidatorDelegate()
{
    delete this->mFactory;
}

QWidget* ValidatorDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QLineEdit *editor = new QLineEdit(parent);
    QValidator* validator = this->mFactory->create(editor);
    editor->setValidator(validator);
    return editor;
}


void ValidatorDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    line->setText(value);
}


void ValidatorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit *line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index, value);
}


void ValidatorDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}
