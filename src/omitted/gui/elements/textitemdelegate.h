#ifndef TEXTITEMDELEGATE_H
#define TEXTITEMDELEGATE_H

#include <QItemDelegate>
#include "textitemwidgetfactory.h"

namespace espreso
{

class TextItemDelegate: public QItemDelegate
{
    Q_OBJECT

public:
    TextItemDelegate(const QString& defaultText,
                     TextItemWidgetFactory* factory,
                     QObject* parent = nullptr);
    ~TextItemDelegate();

protected:
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setEditorData(QWidget * editor, const QModelIndex & index) const override;
    void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
    void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;

private:
    TextItemWidgetFactory* m_factory;
    QString m_defaultText;

private slots:
    void onFinished(QWidget* widget);
};

}

#endif // TEXTITEMDELEGATE_H
