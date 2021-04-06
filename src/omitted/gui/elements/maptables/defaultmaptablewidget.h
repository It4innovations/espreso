#ifndef DEFAULTMAPTABLEWIDGET_H
#define DEFAULTMAPTABLEWIDGET_H

#include "maptablewidget.h"

#include "config/configuration.h"

namespace espreso
{

class DefaultMapTableWidget : public MapTableWidget
{
    Q_OBJECT
public:
    DefaultMapTableWidget(ECFObject* map,
                          const QStringList& headlines,
                          QWidget* parent = 0,
                          std::unique_ptr<CleanRowFactory> rowFactory = std::unique_ptr<CleanRowFactory>(new DefaultCleanRowFactory));

    QVector<QPair<QString, QString> > dataInRows() override;

    bool isValid() override;
    QString errorMessage() override;

    void save() override;

protected:
    ECFObject* m_map;

private:
    void setupColumn(int col, ECFDataType type);
};

}

#endif // DEFAULTMAPTABLEWIDGET_H
