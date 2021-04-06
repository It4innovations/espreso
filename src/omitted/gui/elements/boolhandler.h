#ifndef BOOLHANDLER_H
#define BOOLHANDLER_H

#include <QWidget>

#include "config/configuration.h"

namespace espreso
{

namespace Ui {
class BoolHandler;
}

class BoolHandler : public QWidget
{
    Q_OBJECT

public:
    explicit BoolHandler(ECFParameter*, QWidget *parent = 0);
    ~BoolHandler();

signals:
    void stateChanged();

private slots:
    void onStateChanged(int);

private:
    Ui::BoolHandler *ui;

    ECFParameter* m_param;
};

}

#endif // BOOLHANDLER_H
