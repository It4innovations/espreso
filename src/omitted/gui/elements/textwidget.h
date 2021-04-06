#ifndef TEXTWIDGET_H
#define TEXTWIDGET_H

#include <QWidget>
#include <QLineEdit>

#include "textitemwidget.h"

namespace espreso
{

namespace Ui {
class TextWidget;
}

class TextWidget : public TextItemWidget
{
    Q_OBJECT

public:
    explicit TextWidget(QWidget *parent = 0);
    ~TextWidget();

    void setValidator(QValidator* validator);

    virtual void setText(const QString& text) override;
    virtual QString text() override;

private:
    Ui::TextWidget *ui;
};

}

#endif // TEXTWIDGET_H
