
#ifndef TEXTITEMWIDGET_H
#define TEXTITEMWIDGET_H

#include "iaccessibletext.h"

#include <QWidget>

namespace espreso
{

class TextItemWidget : public QWidget, public IAccessibleText
{
    Q_OBJECT

public:
    explicit TextItemWidget(QWidget* parent = 0) : QWidget(parent) {}
    virtual ~TextItemWidget() {}

    virtual void setText(const QString& text) = 0;
    virtual QString text() = 0;

signals:
    void finished(QWidget* itself);
};

}

#endif // TEXTITEMWIDGET_H
