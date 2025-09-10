
#ifndef OPTIONWIDGET_H
#define OPTIONWIDGET_H

#include "textitemwidget.h"

namespace espreso
{

namespace Ui {
class OptionWidget;
}

class OptionWidget : public TextItemWidget
{
    Q_OBJECT

public:
    explicit OptionWidget(const QStringList& options, QWidget *parent = 0);
    ~OptionWidget();

    void setText(const QString& text) override;
    QString text() override;

private:
    Ui::OptionWidget *ui;

    QStringList m_options;
    QString m_lastOption;

private slots:
    void onIndexChanged(QString);
};

}

#endif // OPTIONWIDGET_H
