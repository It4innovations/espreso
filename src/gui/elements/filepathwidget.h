
#ifndef FILEPATHWIDGET_H
#define FILEPATHWIDGET_H

#include "textitemwidget.h"

namespace espreso
{

namespace Ui {
class FilepathWidget;
}

class FilepathWidget : public TextItemWidget
{
    Q_OBJECT

public:
    explicit FilepathWidget(QWidget *parent = 0);
    ~FilepathWidget();

    virtual void setText(const QString& text) override;
    virtual QString text() override;

private slots:
    void onPressed();

private:
    Ui::FilepathWidget *ui;
};

}

#endif // FILEPATHWIDGET_H
