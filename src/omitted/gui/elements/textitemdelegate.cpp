#include "textitemdelegate.h"

using namespace espreso;

TextItemDelegate::TextItemDelegate(const QString& defaultText,
                                   TextItemWidgetFactory* factory,
                                   QObject *parent)
    : QItemDelegate(parent)
{
    this->m_factory = factory;
    this->m_defaultText = defaultText;
}

TextItemDelegate::~TextItemDelegate()
{
    delete this->m_factory;
}

QWidget* TextItemDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    TextItemWidget *editor = this->m_factory->create(parent);
    editor->setText(this->m_defaultText);
    connect(editor, SIGNAL(finished(QWidget*)), this, SLOT(onFinished(QWidget*)));

    return editor;
}


void TextItemDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    TextItemWidget *widget = static_cast<TextItemWidget*>(editor);
    widget->setText(value);
}


void TextItemDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    TextItemWidget *widget = static_cast<TextItemWidget*>(editor);
    QString value = widget->text();
    model->setData(index, value);
}


void TextItemDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

void TextItemDelegate::onFinished(QWidget *widget)
{
    this->commitData(widget);
}
