#include "ecfparametertreedelegate.h"

#include <QPainter>
#include <QDebug>

#include "textitemwidget.h"

using namespace espreso;

ECFParameterTreeDelegate::ECFParameterTreeDelegate(QObject *parent)
    : QItemDelegate(parent)
{

}

ECFParameterTreeDelegate::~ECFParameterTreeDelegate()
{
    for (auto it = m_editors.begin(); it != m_editors.end(); ++it)
    {
        delete (*it);
    }

    this->m_editors.clear();
}

void ECFParameterTreeDelegate::registerEditor(int editorId, TextItemWidgetFactory *factory)
{
    this->m_editors[editorId] = factory;
}

TextItemWidgetFactory* ECFParameterTreeDelegate::editorFactory(int id)
{
    return this->m_editors[id];
}

void ECFParameterTreeDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    bool ok;
    int editorIndex = index.data(EditorRole).toInt(&ok);
    if (ok
        && this->m_editors[editorIndex] != nullptr
        && dynamic_cast<DataTypeEditWidgetFactory*>(this->m_editors[editorIndex]))
    {
        TextItemWidget *editor = this->m_editors[editorIndex]->create();
		editor->setText(index.data(Qt::EditRole).toString());
        editor->hide();
        editor->resize(option.rect.size());
        QPixmap pixmap(option.rect.size());
        editor->render(&pixmap);
        painter->drawPixmap(option.rect, pixmap);
        delete editor;
    }
    else
    {
        QItemDelegate::paint(painter, option, index);
    }
}

QSize ECFParameterTreeDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    bool ok;
    int editorIndex = index.data(EditorRole).toInt(&ok);
    if (ok
        && this->m_editors[editorIndex] != nullptr
        && dynamic_cast<DataTypeEditWidgetFactory*>(this->m_editors[editorIndex]))
    {
        TextItemWidget *editor = this->m_editors[editorIndex]->create();
        editor->hide();
        QSize hint = editor->sizeHint();
        delete editor;

        return hint;
    }
    else
    {
        return QItemDelegate::sizeHint(option, index);
    }
}

QWidget* ECFParameterTreeDelegate::createEditor(QWidget *parent,
                                    const QStyleOptionViewItem &option,
                                    const QModelIndex &index) const
{
    int editorIndex = index.data(EditorRole).toInt();
    TextItemWidget *editor = this->m_editors[editorIndex]->create(parent);
    connect(editor, SIGNAL(finished(QWidget*)), this, SLOT(onFinished(QWidget*)));

    return editor;
}

void ECFParameterTreeDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    QString value = index.model()->data(index, Qt::EditRole).toString();
    TextItemWidget *widget = static_cast<TextItemWidget*>(editor);
    widget->setText(value);
}


void ECFParameterTreeDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    TextItemWidget *widget = static_cast<TextItemWidget*>(editor);
    QString value = widget->text();
    model->setData(index, value);
}


void ECFParameterTreeDelegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

void ECFParameterTreeDelegate::onFinished(QWidget *widget)
{
    this->commitData(widget);
}
