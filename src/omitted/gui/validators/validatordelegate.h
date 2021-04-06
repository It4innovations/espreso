#ifndef VALIDATORDELEGATE_H
#define VALIDATORDELEGATE_H

#include <QItemDelegate>
#include <QValidator>
#include "validatorfactory.h"

class ValidatorDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    ValidatorDelegate(ValidatorFactory* factory, QObject* parent = nullptr);
    ~ValidatorDelegate();

protected:
    QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override;
    void setEditorData(QWidget * editor, const QModelIndex & index) const override;
    void setModelData(QWidget * editor, QAbstractItemModel * model, const QModelIndex & index) const override;
    void updateEditorGeometry(QWidget * editor, const QStyleOptionViewItem & option, const QModelIndex & index) const override;

private:
    ValidatorFactory* mFactory;
};

#endif // VALIDATORDELEGATE_H
