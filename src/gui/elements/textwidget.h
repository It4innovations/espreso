
#ifndef TEXTWIDGET_H
#define TEXTWIDGET_H

#include "textitemwidget.h"

#include <QWidget>
#include <QLineEdit>

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

    void setText(const QString& text) override;
    QString text() override;

private:
    Ui::TextWidget *ui;
};

}

#endif // TEXTWIDGET_H
