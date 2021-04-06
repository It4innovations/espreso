#ifndef FIELDHANDLER_H
#define FIELDHANDLER_H

#include <QWidget>

#include "ivalidatableobject.h"
#include "isavableobject.h"
#include "validators/validatorfactory.h"
#include "config/configuration.h"

namespace espreso
{

namespace Ui {
class FieldHandler;
}

class FieldHandler : public QWidget, public ISavableObject, public IValidatableObject
{
    Q_OBJECT

public:
    explicit FieldHandler(ECFParameter* data,
                          const ValidatorFactory* validator = nullptr,
                          bool withLabel = false,
                          QWidget *parent = 0);
    ~FieldHandler();

    void setValue(const QString& val);
    QString value() const;

    void save() override;
    void saveState() override;
    void restoreState() override;

    bool isValid() override;
    QString errorMessage() override;

private:
    Ui::FieldHandler *ui;

    ECFParameter* m_data;
    std::string m_saved_state;
};

}

#endif // FIELDHANDLER_H
