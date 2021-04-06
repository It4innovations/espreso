#ifndef PATHHANDLER_H
#define PATHHANDLER_H

#include <QWidget>
#include <QLabel>

#include "config/configuration.h"
#include "ivalidatableobject.h"

namespace espreso
{

class PathHandler : public QWidget, public IValidatableObject
{
    Q_OBJECT

public:
    PathHandler(ECFParameter* path, QWidget* parent = 0);

    bool isValid() override;
    QString errorMessage() override;

private slots:
    void onPressed();

private:
    ECFParameter* m_path;
    QLabel* m_label;
};

}
#endif // PATHHANDLER_H
