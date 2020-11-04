
#ifndef PATHHANDLER_H
#define PATHHANDLER_H

#include "ivalidatableobject.h"

#include "config/configuration.h"

#include <QWidget>
#include <QLabel>

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
