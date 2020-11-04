
#ifndef SPINNERHANDLER_H
#define SPINNERHANDLER_H

#include "ivalidatableobject.h"
#include "isavableobject.h"

#include "config/configuration.h"

#include <QWidget>

namespace espreso
{

namespace Ui {
class SpinnerHandler;
}

class SpinnerHandler : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

public:
    explicit SpinnerHandler(ECFParameter* data,
                            bool withLabel = false,
                            QWidget *parent = 0);
    ~SpinnerHandler();

    void setValue(const QString& val);

    void save() override;
    void saveState() override;
    void restoreState() override;

    bool isValid() override;
    QString errorMessage() override;

signals:
    void valueChanged(int val);

private slots:
    void onSpinnerValueChanged(int val);

private:
    Ui::SpinnerHandler *ui;

    ECFParameter* m_data;
    std::string m_saved_state;
};

}

#endif // SPINNERHANDLER_H
