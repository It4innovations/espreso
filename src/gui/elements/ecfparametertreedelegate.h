
#ifndef ECFPARAMETERTREEDELEGATE_H
#define ECFPARAMETERTREEDELEGATE_H

#include "textitemwidgetfactory.h"

#include <QItemDelegate>
#include <QHash>

namespace espreso
{

class ECFParameterTreeDelegate: public QItemDelegate
{
    Q_OBJECT

public:
    static const int EditorRole = Qt::UserRole + 1;

    ECFParameterTreeDelegate(QObject* parent = nullptr);
    ~ECFParameterTreeDelegate();

    void registerEditor(int editorId, TextItemWidgetFactory *factory);
    TextItemWidgetFactory* editorFactory(int id);

    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    QSize sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const override;

protected:
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setEditorData(QWidget * editor, const QModelIndex & index) const override;
    void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
    void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;

private:
    QString m_defaultText;
    QHash<int, TextItemWidgetFactory*> m_editors;

private slots:
    void onFinished(QWidget* widget);

};

}

#endif // ECFPARAMETERTREEDELEGATE_H
