#ifndef EXPRESSIONEDITDELEGATE_H
#define EXPRESSIONEDITDELEGATE_H

#include "expressionedit.h"

#include <QItemDelegate>
#include <QLineEdit>
#include <QPainter>

namespace espreso
{
    class ExpressionEditDelegate: public QItemDelegate
    {
        Q_OBJECT

    public:
        ExpressionEditDelegate(const std::vector<std::string>& variables, QObject* parent = nullptr);

    signals:
        void validStateChanged(bool valid);

    protected:
        QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
        void setEditorData(QWidget * editor, const QModelIndex & index) const override;
        void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
        void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;
        void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const override;

    private:
        std::vector<std::string> m_variables;
    };
}

#endif // EXPRESSIONEDITDELEGATE_H
