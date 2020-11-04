#ifndef TABLETYPEWIDGET_H
#define TABLETYPEWIDGET_H

#include "tablewidget.h"

namespace espreso
{
    class TableTypeWidget : public TableWidget
    {
        Q_OBJECT

    public:
        static QStringList headlines();

        TableTypeWidget(QWidget *parent = 0);

        virtual void addData(const QString& data) override;

        QString data() override;
    };
}

#endif // TABLETYPEWIDGET_H
