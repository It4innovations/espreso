#include "comboboxdelegate.h"

ComboBoxDelegate::ComboBoxDelegate(const QStringList& options, QObject *parent) : QItemDelegate(parent)
{
    this->mOptions = options;
}

QWidget* ComboBoxDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    QComboBox *editor = new QComboBox(parent);
    editor->addItems(this->mOptions);

    return editor;
}


void ComboBoxDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    QComboBox *cmb = static_cast<QComboBox*>(editor);
    cmb->setCurrentIndex(mOptions.indexOf(value));
}


void ComboBoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QComboBox *cmb = static_cast<QComboBox*>(editor);
    QString value = cmb->currentText();
    model->setData(index, value);
}


void ComboBoxDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}
